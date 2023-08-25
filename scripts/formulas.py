import argparse
import math
from itertools import combinations
import dendropy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.special import lambertw


def is_balanced(st_r, a):
    mrca_123 = st_r.mrca(taxon_labels=[a[0], a[1], a[2]])
    mrca_124 = st_r.mrca(taxon_labels=[a[0], a[1], a[3]])
    mrca_134 = st_r.mrca(taxon_labels=[a[0], a[2], a[3]])
    mrca_234 = st_r.mrca(taxon_labels=[a[1], a[2], a[3]])
    if mrca_123 == mrca_124 == mrca_134 == mrca_234:
        return True
    return False


def get_branch_lengths_rooted(t, labels, balanced, nc=1.0):
    a, b, c, d = labels
    if balanced:
        l_t1 = t.mrca(taxon_labels=[a, b]).edge.length
        l_t2 = t.mrca(taxon_labels=[c, d]).edge.length
    else:
        l_t1 = t.mrca(taxon_labels=[a, b]).edge.length
        l_t2 = t.mrca(taxon_labels=[a, b, c]).edge.length
    l_a = t.find_node_with_taxon_label(a).edge.length
    l_b = t.find_node_with_taxon_label(b).edge.length
    l_c = t.find_node_with_taxon_label(c).edge.length
    l_d = t.find_node_with_taxon_label(d).edge.length
    return l_t1 * nc, l_t2 * nc, l_a * nc, l_b * nc, l_c * nc, l_d * nc


def get_branch_lengths(t, labels, nc=1.0):
    a, b, c, d = labels
    l_a = t.find_node_with_taxon_label(a).edge.length
    l_b = t.find_node_with_taxon_label(b).edge.length
    l_c = t.find_node_with_taxon_label(c).edge.length
    l_d = t.find_node_with_taxon_label(d).edge.length
    l_total = sum(e.length for e in t.postorder_edge_iter() if e.tail_node)
    l_i = l_total - (l_a + l_b + l_c + l_d)
    return l_i * nc, l_a * nc, l_b * nc, l_c * nc, l_d * nc


def node_labels(t):
    a, b, c, d = None, None, None, None
    mrca = [[None for x in range(4)] for y in range(4)]
    labels = ['1', '2', '3', '4']
    for x, y in combinations(labels, 2):
        mrca[labels.index(x)][labels.index(y)] = mrca[labels.index(y)][labels.index(x)] = t.mrca(taxon_labels=[x, y])
    if is_balanced(t, labels):
        a = '1'
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
    print('a,b,c,d:', a, b, c, d)
    return a, b, c, d


def get_mutation_rates(t, labels, rates_path, global_subs_path, balanced):
    # read mutation rates in pre-order traversal of the tree
    for node in t.postorder_node_iter():
        if node.is_leaf():
            node.label = node.taxon.label
        else:
            left = node._child_nodes[0]
            right = node._child_nodes[1]
            node_list = left.label.split(',') + right.label.split(',')
            node.label = ','.join(sorted(node_list))
    preorder_labels = [node.label for node in t.preorder_node_iter()]
    # print(preorder_labels)
    with open(global_subs_path) as f:
        global_rate = float(f.readline().split()[-1])
    print('global subs rate:', global_rate)
    mu_rates = []
    for line in open(rates_path):
        mu_rates = line.split("\t")
        mu_rates = [float(mu) * global_rate for mu in mu_rates[:-1]]
    print('max/min substitution rate:', max(mu_rates) / min(mu_rates))
    if balanced:
        mu1 = mu_rates[preorder_labels.index(','.join(sorted(labels[0:2])))]
        mu2 = mu_rates[preorder_labels.index(','.join(sorted(labels[2:4])))]
        mu3 = mu_rates[preorder_labels.index(','.join(sorted(labels[0:4])))]
        mua = mu_rates[preorder_labels.index(labels[0])]
        mub = mu_rates[preorder_labels.index(labels[1])]
        muc = mu_rates[preorder_labels.index(labels[2])]
        mud = mu_rates[preorder_labels.index(labels[3])]
    else:
        mu1 = mu_rates[preorder_labels.index(','.join(sorted(labels[0:2])))]
        mu2 = mu_rates[preorder_labels.index(','.join(sorted(labels[0:3])))]
        mu3 = mu_rates[preorder_labels.index(','.join(sorted(labels[0:4])))]
        mua = mu_rates[preorder_labels.index(labels[0])]
        mub = mu_rates[preorder_labels.index(labels[1])]
        muc = mu_rates[preorder_labels.index(labels[2])]
        mud = mu_rates[preorder_labels.index(labels[3])]
    print("mu3, mu2, mu1:", mu3, mu2, mu1)
    print("mua, mub, muc, mud:", mua, mub, muc, mud)
    return mu3, mu2, mu1, mua, mub, muc, mud


def balanced_internal_formulas(pop_size, t1, t2, mu1, mu2, mu3):
    lm_i_theory = pop_size * (3 * math.exp(-t1) * (mu1 - mu3) + mu3 * math.exp(-(t1 + t2)) + 3 * math.exp(-t2) * (
            mu2 - mu3) + 3 * (t1 * mu1 + t2 * mu2 - mu1 - mu2 + 2 * mu3)) / (3 - 2 * math.exp(-(t1 + t2)))
    ln_i_theory = pop_size * mu3
    lm_i_simplified = pop_size * ((3 * mu1 * (t1 + math.exp(-t1) - 1)) / (3 - 2 * math.exp(-t1)) + mu3)
    ln_i_simplified = pop_size * mu3
    return lm_i_theory, ln_i_theory, lm_i_simplified, ln_i_simplified


def balanced_terminal_formulas(pop_size, t1, t2, t_a, mu1, mu3, mua):
    # terminal A, B, C, D
    lm_a_theory = pop_size * ((math.exp(-(t1 + t2)) * (-6 * t1 * mu1 - 7 * mu3) + 9 * (
            (1 - math.exp(-t1)) * mu1 + mu3 * math.exp(-t1))) / (9 - 6 * math.exp(-(t1 + t2))) + t_a * mua)
    ln_a_theory = pop_size * (t1 * mu1 + 2 / 3 * mu3 + t_a * mua)
    lm_a_simplified = pop_size * (mu1 + mu1 * (1 + 6 * t1) / (6 - 9 * math.exp(t1 + t2)) + t_a * mua)
    ln_a_simplified = pop_size * (t1 * mu1 + 2 / 3 * mu1 + t_a * mua)
    # delta_a_theory = pop_size * ()
    # delta_a_simplified = pop_size * ()
    return lm_a_theory, ln_a_theory, lm_a_simplified, ln_a_simplified


def unbalanced_internal_formulas(pop_size, t1, t2, mu1, mu2, mu3):
    lm_i_theory = pop_size * (((math.exp(-3 * t2) + 3 * math.exp(-t2) - 6 * math.exp(t1 - t2)) * (mu2 - mu3) + 6 * (
            1 - math.exp(t1) + t1 * math.exp(t1)) * mu1) / (2 * (3 * math.exp(t1) - 2)) + mu2)
    ln_i_theory = pop_size * (mu2 + 0.5 * (mu2 - mu3) * (math.exp(-3 * t2) - 3 * math.exp(-t2)))
    lm_i_simplified = pop_size * ((3 * mu1 * (t1 + math.exp(-t1) - 1)) / (3 - 2 * math.exp(-t1)) + mu2)
    ln_i_simplified = pop_size * mu2
    return lm_i_theory, ln_i_theory, lm_i_simplified, ln_i_simplified


def unbalanced_cherry_formulas(pop_size, t1, t2, t_a, mu1, mu2, mu3, mua):
    # terminal A, B
    lm_a_theory = pop_size * ((6 * t1 * mu1 + 3 * mu1 - mu2 + math.exp(-3 * t2) * (mu2 - 2 * mu3)) / (
            6 - 9 * math.exp(t1)) + mu1 + mua * t_a)
    ln_a_theory = pop_size * (1 / 12 * (10 * mu2 - 9 * math.exp(-t2) * (mu2 - mu3) - 3 * math.exp(-3 * t2) * (
            mu2 + mu3)) + t1 * mu1 + t_a * mua)
    lm_a_simplified = pop_size * ((6 * t1 * mu1 + 3 * mu1 - mu2) / (6 - 9 * math.exp(t1)) + mu1 + t_a * mua)
    ln_a_simplified = pop_size * (t1 * mu1 + 5 / 6 * mu2 + t_a * mua)
    return lm_a_theory, ln_a_theory, lm_a_simplified, ln_a_simplified


def unbalanced_middle_branch_formulas(pop_size, t1, t2, t_c, mu2, mu3, muc):
    # terminal C
    lm_c_theory = pop_size * (-math.exp(-t2) * (mu2 - mu3) + mu2 + muc * t_c + (
            2 * mu2 - (3 * math.exp(-t2) - math.exp(-3 * t2)) * (mu2 - mu3) - 4 * mu3 * math.exp(-3 * t2)) / (
                                      6 * (3 * math.exp(t1) - 2)))
    ln_c_theory = pop_size * (1 / 3 * mu2 * (1 + math.exp(-3 * t2)) + muc * t_c)
    lm_c_simplified = pop_size * (mu2 + muc * t_c + mu2 / (3 * (3 * math.exp(t1) - 2)))
    ln_c_simplified = pop_size * (1 / 3 * mu2 + muc * t_c)
    return lm_c_theory, ln_c_theory, lm_c_simplified, ln_c_simplified


def unbalanced_root_adjacent_formulas(pop_size, t1, t2, t_d, mu2, mu3, mud):
    # terminal D
    lm_d_theory = pop_size * (math.exp(-t2) * (mu2 - mu3) - mu2 + 2 * mu3 + t2 * mu2 + t_d * mud + (
            -2 * mu2 + (3 * math.exp(-t2) - math.exp(-3 * t2)) * (mu2 - mu3)) / (6 * (3 * math.exp(t1) - 2)))
    ln_d_theory = pop_size * ((3 / 2 * math.exp(-t2) - 1 / 6 * math.exp(-3 * t2)) * (
            mu2 - mu3) - 4 / 3 * mu2 + 2 * mu3 + t2 * mu2 + mud * t_d)
    lm_d_simplified = pop_size * (mu2 + t2 * mu2 + mud * t_d - (mu2 / (3 * (3 * math.exp(t1) - 2))))
    ln_d_simplified = pop_size * (2 / 3 * mu2 + t2 * mu2 + mud * t_d)
    return lm_d_theory, ln_d_theory, lm_d_simplified, ln_d_simplified


def compute_su_branches(dataset_name, pop_size, dir, idx, df):
    tns = dendropy.TaxonNamespace()
    st_r = dendropy.Tree.get(path=dir + '/' + idx + '/s_tree.trees', schema='newick', taxon_namespace=tns,
                             rooting="force-rooted")
    st_r_c = dendropy.Tree.get(path=dir + '/' + idx + '/species_trees.cu', schema='newick', taxon_namespace=tns,
                               rooting="force-rooted")
    st_u = dendropy.Tree.get(path=dir + '/' + idx + '/s_tree.trees', schema='newick', taxon_namespace=tns)
    st_c = dendropy.Tree.get(path=dir + '/' + idx + '/species_trees.cu', schema='newick', taxon_namespace=tns)
    gts = dendropy.TreeList.get(path=dir + '/' + idx + '/truegenetrees', schema='newick', taxon_namespace=tns)

    st_c.deroot()
    st_u.deroot()
    for gt in gts:
        gt.deroot()

    print('st_r', st_r)
    print('st_c', st_c)
    print('st_u', st_u)

    balanced = is_balanced(st_r, ['1', '2', '3', '4'])
    if balanced:
        print('balanced')
    else:
        print('unbalanced')
    labels = node_labels(st_r)
    t1, t2, t_a, t_b, t_c, t_d = get_branch_lengths_rooted(st_r_c, labels, balanced, nc=1 / pop_size)
    print('t1, t2, t_a, t_b, t_c, t_d')
    print(t1, t2, t_a, t_b, t_c, t_d)
    mu3, mu2, mu1, mua, mub, muc, mud = get_mutation_rates(st_r, labels, dir + '/' + idx + '/s_tree.ralpha',
                                                           dir + '/' + idx + '/mu.txt', balanced)
    m_gts = []
    n_gts = []
    for gt in gts:
        if dendropy.calculate.treecompare.symmetric_difference(st_u, gt) == 0:
            m_gts.append(gt)
        else:
            n_gts.append(gt)
    num_m_gts = len(m_gts)
    num_n_gts = len(n_gts)
    num_gts = len(gts)
    print('#matching, #non-matching: ', num_m_gts, num_n_gts)
    print('mu3/mu2, mu3/mu1', mu3 / mu2, mu3 / mu1)
    l_i_true, l_a_true, l_b_true, l_c_true, l_d_true = get_branch_lengths(st_u, labels)
    print('true l_a', l_a_true)

    d, _, _, _, _ = get_branch_lengths(st_c, labels, nc=1 / pop_size)
    print('d, t1, t2', d, t1, t2)
    print('ta*mua / l_a_true', t_a * mua * pop_size / l_a_true)
    print('tb*mub / l_b_true', t_b * mub * pop_size / l_b_true)
    print('tc*muc / l_c_true', t_c * muc * pop_size / l_c_true)
    if balanced:
        print('td*mud / l_d_true', t_d * mud * pop_size / l_d_true)
    else:
        print('td*mud+t2mu2 / l_d_true', (t_d * mud + t2 * mu2) * pop_size / l_d_true)
    p = 1 - math.exp(-d)
    print("d, p", d, p)

    p_est = (num_m_gts - 0.5 * (1 + num_n_gts)) / (num_m_gts + num_n_gts + 1)
    d_est = -np.log(1 - p_est)
    print("d_est, p_est", d_est, p_est)

    bl_m_gts = np.zeros((5, len(m_gts)))
    bl_n_gts = np.zeros((5, len(n_gts)))
    for i in range(len(m_gts)):
        bl_m_gts[0][i], bl_m_gts[1][i], bl_m_gts[2][i], bl_m_gts[3][i], bl_m_gts[4][i] = get_branch_lengths(m_gts[i],
                                                                                                            labels)
    for i in range(len(n_gts)):
        bl_n_gts[0][i], bl_n_gts[1][i], bl_n_gts[2][i], bl_n_gts[3][i], bl_n_gts[4][i] = get_branch_lengths(n_gts[i],
                                                                                                            labels)

    l_a_naive = (np.sum(bl_m_gts[1]) + np.sum(bl_n_gts[1])) / num_gts
    l_b_naive = (np.sum(bl_m_gts[2]) + np.sum(bl_n_gts[2])) / num_gts
    l_c_naive = (np.sum(bl_m_gts[3]) + np.sum(bl_n_gts[3])) / num_gts
    l_d_naive = (np.sum(bl_m_gts[4]) + np.sum(bl_n_gts[4])) / num_gts

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

    delta = (lm_i - ln_i) / ln_i if ((lm_i - ln_i) >= 0 and ln_i > 0) else 1e-03
    l_i_est = (delta + np.real(lambertw(-1 / 3 * np.exp(-delta - 1) * (2 * delta + 3))) + 1) * ln_i

    if balanced:
        lm_i_theory, ln_i_theory, lm_i_simplified, ln_i_simplified = balanced_internal_formulas(pop_size, t1, t2, mu1,
                                                                                                mu2, mu3)
        lm_a_theory, ln_a_theory, lm_a_simplified, ln_a_simplified = balanced_terminal_formulas(pop_size, t1, t2, t_a,
                                                                                                mu1, mu3, mua)
        lm_b_theory, ln_b_theory, lm_b_simplified, ln_b_simplified = balanced_terminal_formulas(pop_size, t1, t2, t_b,
                                                                                                mu1, mu3, mub)
        lm_c_theory, ln_c_theory, lm_c_simplified, ln_c_simplified = balanced_terminal_formulas(pop_size, t2, t1, t_c,
                                                                                                mu2, mu3, muc)
        lm_d_theory, ln_d_theory, lm_d_simplified, ln_d_simplified = balanced_terminal_formulas(pop_size, t2, t1, t_d,
                                                                                                mu2, mu3, mud)
        l_a_theory = ln_a_theory - 2 / 3 * mu1 * pop_size - 1 / 3 * (
                mu1 * pop_size * (1 - math.exp(-(t1 + t2))) - (lm_a_theory - ln_a_theory) * (
                3 - 2 * math.exp(-(t1 + t2))))
        l_b_theory = ln_b_theory - 2 / 3 * mu1 * pop_size - 1 / 3 * (
                mu1 * pop_size * (1 - math.exp(-(t1 + t2))) - (lm_b_theory - ln_b_theory) * (
                3 - 2 * math.exp(-(t1 + t2))))
        l_c_theory = ln_c_theory - 2 / 3 * mu2 * pop_size - 1 / 3 * (
                mu2 * pop_size * (1 - math.exp(-(t1 + t2))) - (lm_c_theory - ln_c_theory) * (
                3 - 2 * math.exp(-(t1 + t2))))
        l_d_theory = ln_d_theory - 2 / 3 * mu2 * pop_size - 1 / 3 * (
                mu2 * pop_size * (1 - math.exp(-(t1 + t2))) - (lm_d_theory - ln_d_theory) * (
                3 - 2 * math.exp(-(t1 + t2))))

        mu1_est = ln_i #mu1 * pop_size
        mu2_est = ln_i #mu2 * pop_size
        l_a_est = ln_a - 2 / 3 * mu1_est - 1 / 3 * (mu1_est * p_est - (lm_a - ln_a) * (1 + 2*p_est))
        l_b_est = ln_b - 2 / 3 * mu1_est - 1 / 3 * (mu1_est * p_est - (lm_b - ln_b) * (1 + 2*p_est))
        l_c_est = ln_c - 2 / 3 * mu2_est - 1 / 3 * (mu2_est * p_est - (lm_c - ln_c) * (1 + 2*p_est))
        l_d_est = ln_d - 2 / 3 * mu2_est - 1 / 3 * (mu2_est * p_est - (lm_d - ln_d) * (1 + 2*p_est))

    else:
        lm_i_theory, ln_i_theory, lm_i_simplified, ln_i_simplified = unbalanced_internal_formulas(pop_size, t1, t2, mu1,
                                                                                                  mu2, mu3)
        lm_a_theory, ln_a_theory, lm_a_simplified, ln_a_simplified = unbalanced_cherry_formulas(pop_size, t1, t2, t_a,
                                                                                                mu1, mu2, mu3, mua)
        lm_b_theory, ln_b_theory, lm_b_simplified, ln_b_simplified = unbalanced_cherry_formulas(pop_size, t1, t2, t_b,
                                                                                                mu1, mu2, mu3, mub)
        lm_c_theory, ln_c_theory, lm_c_simplified, ln_c_simplified = unbalanced_middle_branch_formulas(pop_size, t1, t2,
                                                                                                       t_c, mu2, mu3,
                                                                                                       muc)
        lm_d_theory, ln_d_theory, lm_d_simplified, ln_d_simplified = unbalanced_root_adjacent_formulas(pop_size, t1, t2,
                                                                                                       t_d, mu2, mu3,
                                                                                                       mud)
        # the formula based on true values of the parameters as input
        l_a_theory = ln_a_theory + (pop_size * mu1 * (math.exp(-t1) - 1 + t1) + (lm_a_theory - ln_a_theory) * (
                1 - 2 / 3 * math.exp(-t1))) / (1 - 4 / 5 * math.exp(-t1)) - t1 * mu1 * pop_size # should this be l_i_theory?
        l_b_theory = ln_b_theory + (pop_size * mu1 * (math.exp(-t1) - 1 + t1) + (lm_b_theory - ln_b_theory) * (
                    1 - 2 / 3 * math.exp(-t1))) / (1 - 4 / 5 * math.exp(-t1)) - t1 * mu1 * pop_size
        l_c_theory = ln_c_theory - 1 / 3 * (2 - (1 / (2 - math.exp(-t1)))) * (lm_c_theory - ln_c_theory)
        l_d_theory = ln_d_theory - 2 / 3 * (2 + (1 / (1 - math.exp(-t1)))) * (lm_d_theory - ln_d_theory)

        # the formulas based on estimated values
        mu1_est = ln_i #mu1 * pop_size
        l_a_est = ln_a + (mu1_est * (d_est - p_est) + (lm_a - ln_a) * (
                1 - 2 / 3 * (1 - p_est))) / (1 - 4 / 5 * (1 - p_est)) - l_i_est
        l_b_est = ln_b + (mu1_est * (d_est - p_est) + (lm_b - ln_b) * (
                1 - 2 / 3 * (1 - p_est))) / (1 - 4 / 5 * (1 - p_est)) - l_i_est
        l_c_est = ln_c - 1 / 3 * (2 - (1 / (1 + p_est))) * (lm_c - ln_c)
        l_d_est = ln_d - 2 / 3 * (2 + 1 / p_est) * (lm_d - ln_d)

    delta_theory = (lm_i_theory - ln_i_theory) / ln_i_theory if (lm_i_theory - ln_i_theory) >= 0 else 1e-03
    l_i_theory = (delta_theory + np.real(
        lambertw(-1 / 3 * np.exp(-delta_theory - 1) * (2 * delta_theory + 3))) + 1) * ln_i_theory
    l_i_taylor0_1 = ln_i * (1 / 2 * delta + 1 / 6 * math.sqrt(3 * delta * (3 * delta + 4)))
    l_i_est_coal = 1 / 3 * np.abs(lm_i - ln_i) * d_est * (1 + 2 * p_est) / (d_est - p_est)
    l_i_est_lambert = (delta + np.real(lambertw(-1 / 3 * np.exp(-delta - 1) * (2 * delta + 3))) + 1) * ln_i
    l_i_taylor0_2 = (np.sqrt(2 / 3 * delta) + 7 / 9 * delta + 19 * (delta ** (3 / 2)) / (27 * np.sqrt(6)) - (
            delta ** 2) / 9) * ln_i


    print('lm, lm (theory), lm (simplified)', lm_i, lm_i_theory, lm_i_simplified)
    print('ln, ln (theory), ln (simplified)', ln_i, ln_i_theory, ln_i_simplified)
    print('lm(simplified)/lm(theory)', lm_i_simplified / lm_i_theory)
    print('ln(simplified)/ln(theory)', ln_i_simplified / ln_i_theory)
    print('l_i, l_i_est, l_i_theory', l_i_true, l_i_est, l_i_theory)
    print('l_i_theory/l_i', l_i_theory / l_i_true)
    print('l_i_coal/l_i', l_i_est_coal / l_i_true)
    print('l_i_lambert/l_i', l_i_est_lambert / l_i_true)
    print('l_i_taylor0_1/l_i', l_i_taylor0_1 / l_i_true)
    print('l_i_taylor0_2/l_i', l_i_taylor0_2 / l_i_true)
    print('delta, delta_theory, delta/delta_theory', delta, delta_theory, delta/delta_theory)

    print('lm_a, lm_a (theory), lm_a (simplified)', lm_a, lm_a_theory, lm_a_simplified)
    print('ln_a, ln_a (theory), ln_a (simplified)', ln_a, ln_a_theory, ln_a_simplified)
    print('lm_a(simplified)/lm_a(theory)', lm_a_simplified / lm_a_theory)
    print('ln_a(simplified)/ln_a(theory)', ln_a_simplified / ln_a_theory)
    print('l_a_theory/l_a', l_a_theory / l_a_true)

    print('lm_b, lm_b (theory), lm_b (simplified)', lm_b, lm_b_theory, lm_b_simplified)
    print('ln_b, ln_b (theory), ln_b (simplified)', ln_b, ln_b_theory, ln_b_simplified)
    print('lm_b(simplified)/lm_b(theory)', lm_b_simplified / lm_b_theory)
    print('ln_b(simplified)/ln_b(theory)', ln_b_simplified / ln_b_theory)
    print('l_b_theory/l_b', l_b_theory / l_b_true)

    print('lm_c, lm_c (theory), lm_c (simplified)', lm_c, lm_c_theory, lm_c_simplified)
    print('ln_c, ln_c (theory), ln_c (simplified)', ln_c, ln_c_theory, ln_c_simplified)
    print('lm_c(simplified)/lm_c(theory)', lm_c_simplified / lm_c_theory)
    print('ln_c(simplified)/ln_c(theory)', ln_c_simplified / ln_c_theory)
    print('l_c_theory/l_c', l_c_theory / l_c_true)

    print('lm_d, lm_d (theory), lm_d (simplified)', lm_d, lm_d_theory, lm_d_simplified)
    print('ln_d, ln_d (thoery), ln_d (simplified)', ln_d, ln_d_theory, ln_d_simplified)
    print('lm_d(simplified)/lm_d(theory)', lm_d_simplified / lm_d_theory)
    print('ln_d(simplified)/ln_d(theory)', ln_d_simplified / ln_d_theory)
    print('l_d_theory/l_d', l_d_theory / l_d_true)

    df.loc[len(df.index)] = [dataset_name, 'internal', 'balanced' if balanced else 'unbalanced',
                             t1, t2, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             lm_i_theory, lm_i_simplified, lm_i, 'matching',
                             l_i_true, l_i_theory, l_i_est]#, l_i_est_lambert, l_i_est_coal, l_i_taylor0_1, l_i_taylor0_2]

    df.loc[len(df.index)] = [dataset_name, 'internal', 'balanced' if balanced else 'unbalanced',
                             t1, t2, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             ln_i_theory, ln_i_simplified, ln_i, 'non-matching',
                             l_i_true, l_i_theory, l_i_est]#, l_i_est_lambert, l_i_est_coal, l_i_taylor0_1, l_i_taylor0_2]

    df.loc[len(df.index)] = [dataset_name, 'terminal A', 'balanced' if balanced else 'unbalanced',
                             t1, t2, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             lm_a_theory, lm_a_simplified, lm_a, 'matching',
                             l_a_true, l_a_theory, l_a_est]#, 0, 0, 0, 0]

    df.loc[len(df.index)] = [dataset_name, 'terminal A', 'balanced' if balanced else 'unbalanced',
                             t1, t2, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             ln_a_theory, ln_a_simplified, ln_a, 'non-matching',
                             l_a_true, l_a_theory, l_a_est]#, 0, 0, 0, 0]

    df.loc[len(df.index)] = [dataset_name, 'terminal B', 'balanced' if balanced else 'unbalanced',
                             t1, t2, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             lm_b_theory, lm_b_simplified, lm_b, 'matching',
                             l_b_true, l_b_theory, l_b_est]#, 0, 0, 0, 0]

    df.loc[len(df.index)] = [dataset_name, 'terminal B', 'balanced' if balanced else 'unbalanced',
                             t1, t2, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             ln_b_theory, ln_b_simplified, ln_b, 'non-matching',
                             l_b_true, l_b_theory, l_b_est]#, 0, 0, 0, 0]

    df.loc[len(df.index)] = [dataset_name, 'terminal C', 'balanced' if balanced else 'unbalanced',
                             t1, t2, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             lm_c_theory, lm_c_simplified, lm_c, 'matching',
                             l_c_true, l_c_theory, l_c_est]#, 0, 0, 0, 0]

    df.loc[len(df.index)] = [dataset_name, 'terminal C', 'balanced' if balanced else 'unbalanced',
                             t1, t2, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             ln_c_theory, ln_c_simplified, ln_c, 'non-matching',
                             l_c_true, l_c_theory, l_c_est]#, 0, 0, 0, 0]

    df.loc[len(df.index)] = [dataset_name, 'terminal D', 'balanced' if balanced else 'unbalanced',
                             t1, t2, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             lm_d_theory, lm_d_simplified, lm_d, 'matching',
                             l_d_true, l_d_theory, l_d_est]#, 0, 0, 0, 0]

    df.loc[len(df.index)] = [dataset_name, 'terminal D', 'balanced' if balanced else 'unbalanced',
                             t1, t2, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             ln_d_theory, ln_d_simplified, ln_d, 'non-matching',
                             l_d_true, l_d_theory, l_d_est]#, 0, 0, 0, 0]


    '''df.loc[len(df.index)] = ['internal', 'balanced' if balanced else 'unbalanced',
                             d, d_est, p, p_est, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             l_i_true, l_i_lambert_taylor0, l_i_lambert_taylor0/l_i_true,
                             lm_i_theory, lm_i_simplified, lm_i,
                             ln_i_theory, ln_i_simplified , ln_i]

    df.loc[len(df.index)] = ['terminal A', 'balanced' if balanced else 'unbalanced',
                             d, d_est, p, p_est, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             l_a_true, l_a_est, l_a_est/l_a_true,
                             lm_a_theory, lm_a_simplified, lm_a,
                             ln_a_theory, ln_a_simplified, ln_a]

    df.loc[len(df.index)] = ['terminal B', 'balanced' if balanced else 'unbalanced',
                             d, d_est, p, p_est, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             l_b_true, l_b_est, l_b_est/l_b_true,
                             lm_b_theory, lm_b_simplified, lm_b,
                             ln_b_theory, ln_b_simplified, ln_b]

    df.loc[len(df.index)] = ['terminal C', 'balanced' if balanced else 'unbalanced',
                             d, d_est, p, p_est, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             l_c_true, l_c_est, l_c_est/l_c_true,
                             lm_c_theory, lm_c_simplified, lm_c,
                             ln_c_theory, ln_c_simplified, ln_c]

    df.loc[len(df.index)] = ['terminal D', 'balanced' if balanced else 'unbalanced',
                             d, d_est, p, p_est, num_m_gts, num_n_gts,
                             mu3, mu2, mu1, mua, mub, muc, mud,
                             l_d_true, l_d_est, l_d_est/l_d_true,
                             lm_d_theory, lm_d_simplified, lm_d,
                             ln_d_theory, ln_d_simplified, ln_d]'''
    return df


def plot_x_y_line(ax):
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]),
            np.max([ax.get_xlim(), ax.get_ylim()])]
    ax.plot(lims, lims, '--', alpha=0.75, zorder=0, color='black')
    ax.set_aspect('equal')


def plot_correlations(df, name1, name2):
    plt.cla()

    plt.grid(linestyle='--', linewidth=0.5)
    ax = plt.gca()
    ax.set_xlabel('log10 ('+name1+' length)')
    ax.set_ylabel('log10 ('+name2+' length)')
    plot_x_y_line(ax)
    plt.savefig(name1+'_'+name2+'_correlations.pdf', bbox_inches='tight')


def plot_correlations(df, output_name):
    plt.cla()
    #fig, axes = plt.subplots(5, 3, figsize=(18, 24))
    df = df.replace([np.inf, -np.inf], np.nan)
    branches = ['internal', 'terminal A', 'terminal B', 'terminal C', 'terminal D']

    sns.lmplot(data=df, x="l_true", y="l_est", hue="type", row='branch', palette='Dark2',
                    scatter_kws={'alpha': 0.8, 'linewidth': 0}, order=0.5, ci=None)
    #sns.lmplot(data=df, x="lm_theory", y="lm_simplified", hue="type", 'pallete'='Dark2',
    #                scatter_kws={'alpha': 0.8, 'linewidth': 0}, order=2, ci=0)
    #sns.lmplot(data=df, x="ln_theory", y="ln_simplified", hue="type", 'pallete'='Dark2',
    #                scatter_kws={'alpha': 0.8, 'linewidth': 0}, order=2, ci=0)

    axes = plt.gca()
    axes.set_yscale('log')
    axes.set_xscale('log')
    plot_x_y_line(axes)
    '''for ax in axes:
        ax.set_yscale('log')
        ax.set_xscale('log')'''
    ''' ax.set_xlabel('log10 (' + name1 + ' length)')
    ax.set_ylabel('log10 (' + name2 + ' length)')
    ax.text(-0.3, -0.5, 'y=x')
    plot_x_y_line(ax)'''
        
    '''   for j in range(3):
            axes[i, j].text(-0.3, -0.5, 'y=x')
            #axes[i, j].set_title(branches[i])
            # plot_x_y_line(axes[i, j])
            axes[i, j].set_yscale('log')
            axes[i, j].set_xscale('log')'''
        #sns.scatterplot( data=df[df['branch'] == branches[i]], x="l_true", y="l_est", hue="type", palette='Dark2')
        #sns.scatterplot(ax=axes[i, 1], data=df[df['branch'] == branches[i]], x="lm_theory", y="lm_simplified", hue="type", palette='Dark2')
        #sns.scatterplot(ax=axes[i, 2], data=df[df['branch'] == branches[i]], x="ln_theory", y="ln_simplified", hue="type", palette='Dark2')

    plt.savefig(output_name + '.pdf', bbox_inches='tight')


def plot_length_accuracy(df):
    sns.set_theme()
    fig = plt.figure(figsize=(6, 12), dpi=80)
    g = sns.FacetGrid(df, col="branch", hue="type")
    g.map_dataframe(sns.scatterplot, x="log10(l_true)", y="log10(l_est)", palette='viridis')
    g.set_titles(row_template='{row_name}', col_template='{col_name}')
    for ax in g.axes[0]:
        plot_x_y_line(ax)
    g.add_legend()
    plt.savefig(output_name + '_accuracy.pdf', bbox_inches='tight')


def add_additional_columns(df):
    df['log10(mu3/mu2)'] = np.log10(df['mu3/mu2'])
    df['log10(mu3/mu1)'] = np.log10(df['mu3/mu1'])
    df['log10(ln_est/ln_theory)'] = np.log10(df['ln_est/ln_theory'])
    df['log10(lm_est/lm_theory)'] = np.log10(df['lm_est/lm_theory'])
    df['log10(l_est/l_true)'] = np.log10(df['l_est/l_true'])
    df['rmse'] = np.abs(df['l_est'] - df['l_true'])
    df['log-error'] = np.abs(df['log10(l_est)'] - df['log10(l_true)'])
    df['log10(d)'] = np.log10(df['d'])
    df['-log10(d)+|log10(m3/m2)|'] = -np.log10(df['d']) + np.abs(np.log10(df['mu3/mu2']))
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Quartet branch length formulas in substitution units")
    parser.add_argument("-ps", "--popsize", type=str, required=True,
                        help="Population size")
    parser.add_argument("-d", "--datadir", type=str, required=True,
                        help="Quartet data directory")
    parser.add_argument("-n", "--name", type=str, required=True,
                        help="Output file name")
    args = parser.parse_args()
    output_name = args.name + '_formulas'
    df = pd.DataFrame(
        columns=["dataset", "branch", "quartet_type", "t1", "t2", "m_gts", "n_gts",
                 "mu3", "mu2", "mu1", "mua", "mub", "muc", "mud",
                 "l_theory", "l_simplified", "l_empirical", "gene_type",
                 "bl_true", "bl_formula", "bl_est"])#, 'bl_lambert', 'bl_coal', 'bl_taylor0_1', 'bl_taylor0_2'])
    for i in [20, 24, 26, 66, 99, 116, 118, 146]:#range(1, 201):
        print('\n', i)
        try:
            df = compute_su_branches(args.name, int(args.popsize), args.datadir, str(i).zfill(3), df)
        except Exception as e:
            print(str(e))
            continue
    df.to_csv(output_name + '.csv')
    #df = pd.read_csv(output_name + '.csv')
    # df = add_additional_columns(df)
    # plot_length_accuracy(df)
    #plot_correlations(df, output_name)
