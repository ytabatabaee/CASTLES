import argparse
import matplotlib.pyplot as plt
import dendropy
import math
import pandas as pd
import seaborn as sns
import numpy as np


def compare_bl(args):
    tns = dendropy.TaxonNamespace()
    t1 = dendropy.Tree.get(path=args.tree1, schema='newick', taxon_namespace=tns)
    t2 = dendropy.Tree.get(path=args.tree2, schema='newick', taxon_namespace=tns)
    t1.deroot()
    t2.deroot()
    length_diffs = dendropy.calculate.treecompare._get_length_diffs(t1, t2)

    sns.set_theme()
    output_name = 'correlations_'+args.tree2+'.pdf'
    df = pd.DataFrame(columns=['type', "truel", "estl", 'log(truel)', 'log(estl)'])

    squared_error = 0
    log_error = 0
    neg_branches = 0

    idx = 0
    for node in t1.postorder_node_iter():
        (l1, l2) = length_diffs[idx]
        node_type = 'terminal' if node.is_leaf()  else 'internal'
        node_label = node.taxon.label if node.is_leaf()  else ''
        print(node_label, node_type, l1, l2)
        squared_error += (l1 - l2) ** 2
        idx += 1
        if l2 > 0 and l1 > 0:
            log_error += np.abs(math.log10(l1/l2))
            df.loc[len(df.index)] = [node_type, l1, l2, math.log10(l1) if l1 > 0 else 0, math.log10(l2) if l2 > 0 else 0]
        elif l2 < 0:
            neg_branches += 1

    rmse = np.sqrt(squared_error / len(length_diffs))
    log_error /= (len(length_diffs) - neg_branches)
    print('rmse, log_error:', rmse, log_error)
    print('number of negative branches:', neg_branches)

    sns.scatterplot(data=df, x='log(truel)', y='log(estl)', palette='viridis', alpha=0.7, linewidth=0, hue='type')
    #sns.regplot(data=df, x='log(truel)', y='log(estl)', order=2)
    plt.savefig(output_name, bbox_inches='tight')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="compare tree branch lengths")
    parser.add_argument("-t1", "--tree1", type=str,  required=True,
                        help="tree file with branch lengths in newick format")
    parser.add_argument("-t2", "--tree2", type=str, required=True,
                        help="tree file with branch lengths in newick format")
    compare_bl(parser.parse_args())
