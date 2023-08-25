import dendropy
import argparse
import re
import numpy as np

def average_terminal_bl(gts, taxon_label):
    sum_bl = 0
    for gt in gts:
        sum_bl += gt.find_node_with_taxon_label(taxon_label).edge.length
    return sum_bl/len(gts)


def process_node_annotation(annot_str):
    annotations = list(filter(None, re.split(r'[,{}:]', str(annot_str))))
    label_dict = dict()
    label_dict['support'] = float(annotations[1])
    label_dict['length'] = float(annotations[3])
    label_dict['LR_SO'] = dict()
    label_dict['LR_SO']['quartetCnt'] = float(annotations[6])
    label_dict['LR_SO']['sumInternal'] = float(annotations[8])
    label_dict['LR_SO']['sumL'] = float(annotations[10])
    label_dict['LR_SO']['sumR'] = float(annotations[12])
    label_dict['LR_SO']['sumS'] = float(annotations[14])
    label_dict['LR_SO']['sumO'] = float(annotations[16])
    label_dict['LS_RO'] = dict()
    label_dict['LS_RO']['quartetCnt'] = float(annotations[19])
    label_dict['LS_RO']['sumInternal'] = float(annotations[21])
    label_dict['LS_RO']['sumL'] = float(annotations[23])
    label_dict['LS_RO']['sumR'] = float(annotations[25])
    label_dict['LS_RO']['sumS'] = float(annotations[27])
    label_dict['LS_RO']['sumO'] = float(annotations[29])
    label_dict['LO_RS'] = dict()
    label_dict['LO_RS']['quartetCnt'] = float(annotations[32])
    label_dict['LO_RS']['sumInternal'] = float(annotations[34])
    label_dict['LO_RS']['sumL'] = float(annotations[36])
    label_dict['LO_RS']['sumR'] = float(annotations[38])
    label_dict['LO_RS']['sumS'] = float(annotations[40])
    label_dict['LO_RS']['sumO'] = float(annotations[42])
    return label_dict

def safe_div(n, d):
    return n / d if d else 0


def set_branch_length(edge, length):
    edge.length = np.abs(length)
    return edge


def find_sibling(node):
    if node.parent_node._child_nodes[0] == node:
        sibling = node.parent_node._child_nodes[1]
    if node.parent_node._child_nodes[1] == node:
        sibling = node.parent_node._child_nodes[0]
    return sibling


def castles(args):
    # this code is for handling trees with more than 4 taxa, see castles_quartets.py for trees with 4 taxa
    tns = dendropy.TaxonNamespace()
    st = dendropy.Tree.get(path=args.speciestree, schema='newick', taxon_namespace=tns)
    gts = dendropy.TreeList.get(path=args.genetrees, schema='newick', taxon_namespace=tns)

    print('* CASTLES: Coalescent-Aware Species Tree Length Estimation in Substitution-units *\n')

    print('Number of species:', len(tns))
    print('Number of genes:', len(gts))

    print('\nCalculating lengths in post-order traversal of internal nodes:')

    for node in st.postorder_node_iter():
        if node.taxon is not None:
            continue
        else:
            label_dict = process_node_annotation(node.label)
            num_m_gts = label_dict['LR_SO']['quartetCnt']
            num_n_gts = label_dict['LS_RO']['quartetCnt'] + label_dict['LO_RS']['quartetCnt']
            # calculate coalescent unit lenght
            p_est = (num_m_gts - 0.5 * (1 + num_n_gts)) / (num_n_gts + num_m_gts + 1)
            d_est = -np.log(1 - p_est)

            left = node._child_nodes[0]
            right = node._child_nodes[1]
            sibling = None
            if node.parent_node is not None:
                sibling = find_sibling(node)
            else:
                node.label = None
                continue

            # calculating average branch lenghts in matching and non-matching gene trees
            lm_i = safe_div(label_dict['LR_SO']['sumInternal'], label_dict['LR_SO']['quartetCnt'])
            ln_i = safe_div(label_dict['LS_RO']['sumInternal']+label_dict['LO_RS']['sumInternal'],
                            label_dict['LS_RO']['quartetCnt']+ label_dict['LO_RS']['quartetCnt'])
            lm_a = safe_div(label_dict['LR_SO']['sumL'], label_dict['LR_SO']['quartetCnt'])
            ln_a = safe_div(label_dict['LS_RO']['sumL']+label_dict['LO_RS']['sumL'],
                            label_dict['LS_RO']['quartetCnt']+label_dict['LO_RS']['quartetCnt'])
            lm_b = safe_div(label_dict['LR_SO']['sumR'], label_dict['LR_SO']['quartetCnt'])
            ln_b = safe_div(label_dict['LS_RO']['sumR']+label_dict['LO_RS']['sumR'],
                            label_dict['LS_RO']['quartetCnt']+label_dict['LO_RS']['quartetCnt'])
            lm_c = safe_div(label_dict['LR_SO']['sumS'], label_dict['LR_SO']['quartetCnt'])
            ln_c = safe_div(label_dict['LS_RO']['sumS']+label_dict['LO_RS']['sumS'],
                            label_dict['LS_RO']['quartetCnt']+label_dict['LO_RS']['quartetCnt'])
            lm_d = safe_div(label_dict['LR_SO']['sumO'], label_dict['LR_SO']['quartetCnt'])
            ln_d = safe_div(label_dict['LS_RO']['sumO']+label_dict['LO_RS']['sumO'],
                            label_dict['LS_RO']['quartetCnt']+label_dict['LO_RS']['quartetCnt'])

            # compute internal branch (Table S3 in paper)
            delta = safe_div(lm_i - ln_i, ln_i) if lm_i > ln_i else 1e-03
            l_est = 1/6 * (3 * delta + np.sqrt(3 * delta * (4 + 3 * delta))) * ln_i
            mu1_est = ln_i #l_est / d_est

            # compute terminal branches (Table S3 in paper, unbalanced)
            l_a_est = ln_a + (mu1_est * (d_est - p_est) + (lm_a - ln_a) * (
                    1 - 2 / 3 * (1 - p_est))) / (1 - 4 / 5 * (1 - p_est)) - l_est
            l_b_est = ln_b + (mu1_est * (d_est - p_est) + (lm_b - ln_b) * (
                    1 - 2 / 3 * (1 - p_est))) / (1 - 4 / 5 * (1 - p_est)) - l_est
            l_c_est = ln_c - 1 / 3 * (2 - 1 / (p_est + 1)) * (lm_c - ln_c)
            l_d_est = ln_d - 2 / 3 * (2 + 1 / p_est) * (lm_d - ln_d)

            if node.parent_node.parent_node is None: # node is child of the root
                if sibling.is_leaf():
                    if l_d_est == 0:
                        l_d_est = average_terminal_bl(gts, sibling.taxon.label)
                    set_branch_length(node.edge, l_d_est)
                    print('outgroup terminal branch ', sibling.taxon.label, '| length =', l_d_est)
                else:
                    set_branch_length(node.edge, l_est)
                    print('internal root branch | length =', l_est)
            elif not node.edge.length:
                set_branch_length(node.edge, l_est)
                print('internal branch | length =', l_est)
            if left.is_leaf() and not left.edge.length:
                set_branch_length(left.edge, l_a_est)
                print('cherry terminal branch', left.taxon.label,  '| length =', l_a_est)
            if right.is_leaf() and not right.edge.length:
                set_branch_length(right.edge, l_b_est)
                print('cherry terminal branch', right.taxon.label, '| length =', l_b_est)
            if sibling and sibling.is_leaf() and not sibling.edge.length:
                set_branch_length(sibling.edge, l_c_est)
                print('middle terminal branch', sibling.taxon.label, '| length =', l_c_est)
        node.label = None

    st.deroot()
    with open(args.outputtree, 'w') as f:
        f.write(str(st) + ';\n')
        print('\nSpecies tree with SU lengths written to', args.outputtree)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CASTLES")
    parser.add_argument("-t", "--speciestree", type=str,  required=True,
                        help="ASTER-Annotated species tree file in newick format")
    parser.add_argument("-g", "--genetrees", type=str, required=True,
                        help="Gene tree file in newick format")
    parser.add_argument("-o", "--outputtree", type=str, required=True,
                        help="Species tree with SU branch lengths in newick format")
    castles(parser.parse_args())
