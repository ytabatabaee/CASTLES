import dendropy
import argparse
import re
import numpy as np


def average_terminal_bl(gts, taxon_label):
    sum_bl = 0
    for gt in gts:
        sum_bl += gt.find_node_with_taxon_label(taxon_label).edge.length
    return sum_bl/len(gts)


def safe_div(n, d):
    return n / d if d else 0


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


def naive(args):
    tns = dendropy.TaxonNamespace()
    st = dendropy.Tree.get(path=args.speciestree, schema='newick', taxon_namespace=tns)
    gts = dendropy.TreeList.get(path=args.genetrees, schema='newick', taxon_namespace=tns)

    st.deroot()
    for gt in gts:
        gt.deroot()

    for node in st.postorder_node_iter():
        if node.taxon is not None:
            print(node.taxon)
            node.edge.length = average_terminal_bl(gts, node.taxon.label)
        else:
            label_dict = process_node_annotation(node.label)
            num_m_gts = label_dict['LR_SO']['quartetCnt']
            num_n_gts = label_dict['LS_RO']['quartetCnt'] + label_dict['LO_RS']['quartetCnt']
            p_est = (num_m_gts - 0.5 * (num_n_gts+1)) / (num_n_gts + num_m_gts+1)
            d_est = -np.log(1 - p_est)
            print("d_est, p_est", d_est, p_est)
            ln_i = safe_div(label_dict['LS_RO']['sumInternal']+label_dict['LO_RS']['sumInternal'],
                            label_dict['LS_RO']['quartetCnt']+ label_dict['LO_RS']['quartetCnt'])
            node.edge.length = ln_i * d_est
            node.label = None

    with open(args.outputtree, 'w') as f:
        f.write(str(st) + ';\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Naive SU branch length")
    parser.add_argument("-t", "--speciestree", type=str,  required=True,
                        help="Species tree file in newick format")
    parser.add_argument("-g", "--genetrees", type=str, required=True,
                        help="Gene trees file in newick format")
    parser.add_argument("-o", "--outputtree", type=str, required=True,
                        help="Species tree annotated with SU branches in newick format")
    naive(parser.parse_args())
