import dendropy
import argparse
import re
import math
import numpy as np
from itertools import combinations

def is_balanced(st_r, a):
    mrca_123 = st_r.mrca(taxon_labels=[a[0], a[1], a[2]])
    mrca_124 = st_r.mrca(taxon_labels=[a[0], a[1], a[3]])
    mrca_134 = st_r.mrca(taxon_labels=[a[0], a[2], a[3]])
    mrca_234 = st_r.mrca(taxon_labels=[a[1], a[2], a[3]])
    if mrca_123 == mrca_124 == mrca_134 == mrca_234:
        return True
    return False

def get_branch_lengths(t, labels):
    a, b, c, d = labels
    l_a = t.find_node_with_taxon_label(a).edge.length
    l_b = t.find_node_with_taxon_label(b).edge.length
    l_c = t.find_node_with_taxon_label(c).edge.length
    l_d = t.find_node_with_taxon_label(d).edge.length
    l_total = sum(e.length for e in t.postorder_edge_iter() if e.tail_node)
    l_i = l_total - (l_a + l_b + l_c + l_d)
    return l_i, l_a, l_b, l_c, l_d

def node_labels(t, labels):
    a, b, c, d = None, None, None, None
    mrca = [[None for x in range(4)] for y in range(4)]
    for x, y in combinations(labels, 2):
        mrca[labels.index(x)][labels.index(y)] = mrca[labels.index(y)][labels.index(x)] = t.mrca(taxon_labels=[x, y])
    if is_balanced(t, labels):
        a = labels[0]
        for x in range(1, 4):
            if set(mrca[x]) == set(mrca[0]):
                b = labels[x]
            elif not c:
                c = labels[x]
            else:
                d = labels[x]
    else:
        for x in range(4):
            num_unique_mrcas = len(set(mrca[x]))
            if num_unique_mrcas == 2:
                d = labels[x]
            elif num_unique_mrcas == 3:
                c = labels[x]
            elif num_unique_mrcas == 4 and not a:
                a = labels[x]
            else:
                b = labels[x]
    return a, b, c, d


def castles_quartets(st, gts, tns):
    for gt in gts:
        gt.deroot()

    st_u = dendropy.Tree(st)
    st_u.deroot()

    balanced = is_balanced(st, tns.labels())
    labels = node_labels(st, tns.labels())

    m_gts = []
    n_gts = []
    for gt in gts:
        if dendropy.calculate.treecompare.symmetric_difference(st_u, gt) == 0:
            m_gts.append(gt)
        else:
            n_gts.append(gt)
    num_m_gts = len(m_gts)
    num_n_gts = len(n_gts)
    p_est = (num_m_gts - 0.5 * (1 + num_n_gts)) / (num_n_gts + num_m_gts + 1)
    d_est = -np.log(1 - p_est)

    bl_m_gts = np.zeros((5, len(m_gts)))
    bl_n_gts = np.zeros((5, len(n_gts)))
    for i in range(len(m_gts)):
        bl_m_gts[0][i], bl_m_gts[1][i], bl_m_gts[2][i], bl_m_gts[3][i], bl_m_gts[4][i] = get_branch_lengths(m_gts[i], labels)
    for i in range(len(n_gts)):
        bl_n_gts[0][i], bl_n_gts[1][i], bl_n_gts[2][i], bl_n_gts[3][i], bl_n_gts[4][i] = get_branch_lengths(n_gts[i], labels)

    lm_i = np.mean(bl_m_gts[0]) if len(m_gts) > 0 else 0
    ln_i = np.mean(bl_n_gts[0]) if len(n_gts) > 0 else 0
    lm_a = np.mean(bl_m_gts[1]) if len(m_gts) > 0 else 0
    ln_a = np.mean(bl_n_gts[1]) if len(n_gts) > 0 else 0
    lm_b = np.mean(bl_m_gts[2]) if len(m_gts) > 0 else 0
    ln_b = np.mean(bl_n_gts[2]) if len(n_gts) > 0 else 0
    lm_c = np.mean(bl_m_gts[3]) if len(m_gts) > 0 else 0
    ln_c = np.mean(bl_n_gts[3]) if len(n_gts) > 0 else 0
    lm_d = np.mean(bl_m_gts[4]) if len(m_gts) > 0 else 0
    ln_d = np.mean(bl_n_gts[4]) if len(n_gts) > 0 else 0

    delta = (lm_i - ln_i) / ln_i if lm_i > ln_i else 1e-03
    l_est = 1/6 * (3 * delta + np.sqrt(3 * delta * (4 + 3 * delta))) * ln_i
    mu1_est = l_est / d_est

    mu2_est_a = -(mu1_est * 3 * (d_est - p_est) + (lm_a - ln_a) * (1 + 2 * p_est)) * 2 / (1 + 4 * p_est)
    mu2_est_b = -(mu1_est * 3 * (d_est - p_est) + (lm_b - ln_b) * (1 + 2 * p_est)) * 2 / (1 + 4 * p_est)
    mu2_est_c = -(mu1_est * 3 * (d_est - p_est) + (lm_c - ln_c) * (1 + 2 * p_est)) * 2 / (1 + 4 * p_est)
    mu2_est_d = -(mu1_est * 3 * (d_est - p_est) + (lm_d - ln_d) * (1 + 2 * p_est)) * 2 / (1 + 4 * p_est)

    if balanced:
        l_a_est = ln_a - 2 / 3 * mu1_est - 1 / 3 *(mu1_est * p_est - (lm_a - ln_a) * (1 + 2*p_est))
        l_b_est = ln_b - 2 / 3 * mu1_est - 1 / 3 *(mu1_est * p_est - (lm_b - ln_b) * (1 + 2*p_est))
        l_c_est = ln_c - 2 / 3 * mu1_est - 1 / 3 *(mu1_est * p_est - (lm_c - ln_c) * (1 + 2*p_est))
        l_d_est = ln_d - 2 / 3 * mu1_est - 1 / 3 *(mu1_est * p_est - (lm_d - ln_d) * (1 + 2*p_est))
    else:
        l_a_est = ln_a - 5 / 6 * mu2_est_a - l_est
        l_b_est = ln_b - 5 / 6 * mu2_est_b - l_est
        l_c_est = ln_c - 1 / 3 * (2 - 1 / (p_est + 1)) * (lm_c - ln_c)
        l_d_est = ln_d - 2 / 3 * (2 + 1 / p_est) * (lm_d - ln_d)

    st.deroot()

    for node in st.postorder_node_iter():
        if not node.parent_node:
            continue
        if node.taxon is not None:
            if node.taxon.label == labels[0]:
                node.edge.length = l_a_est
                print('cherry terminal branch', node.taxon.label, '| length =', l_a_est)
            elif node.taxon.label == labels[1]:
                node.edge.length = l_b_est
                print('cherry terminal branch', node.taxon.label, '| length =', l_b_est)
            elif node.taxon.label == labels[2]:
                node.edge.length = l_c_est
                if balanced:
                    print('cherry terminal branch', node.taxon.label, '| length =', l_c_est)
                else:
                    print('middle terminal branch', node.taxon.label, '| length =', l_c_est)
            elif node.taxon.label == labels[3]:
                node.edge.length = l_d_est
                if balanced:
                    print('cherry terminal branch', node.taxon.label, '| length =', l_d_est)
                else:
                    print('outgroup terminal branch', node.taxon.label, '| length =', l_d_est)
        else:
            node.edge.length = l_est
        node.label = None

    with open(args.outputtree, 'w') as f:
        f.write(str(st) + ';\n')
        print('\nSpecies tree with SU lengths written to', args.outputtree)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CASTLES")
    parser.add_argument("-t", "--speciestree", type=str,  required=True,
                        help="Species tree in newick format")
    parser.add_argument("-g", "--genetrees", type=str, required=True,
                        help="Gene trees in newick format")
    parser.add_argument("-o", "--outputtree", type=str, required=True,
                        help="Species tree with SU branch lengths in newick format")
    args = parser.parse_args()

    print('* CASTLES: Coalescent-Aware Species Tree Length Estimation in Substitution-units *\n')

    tns = dendropy.TaxonNamespace()
    st = dendropy.Tree.get(path=args.speciestree, schema='newick', taxon_namespace=tns, rooting='force-rooted')
    gts = dendropy.TreeList.get(path=args.genetrees, schema='newick', taxon_namespace=tns)
    castles_quartets(st, gts, tns)