import dendropy
import sys
import numpy as np
import argparse


def main(args):
    tns = dendropy.TaxonNamespace()
    st = dendropy.Tree.get(path=args.speciestree, schema='newick', taxon_namespace=tns)
    gts = dendropy.TreeList.get(path=args.genetrees, schema='newick', taxon_namespace=tns)

    dist_mat = np.zeros((len(tns), len(tns), len(gts)))
    idx2 = 0
    for idx in range(len(gts)):
        gts[idx].deroot()
        pdc = gts[idx].phylogenetic_distance_matrix()
        exclude = True
        for i in range(len(tns)):
            for j in range(i + 1, len(tns)):
                dist_mat[i][j][idx2] = pdc(tns[i], tns[j])
                dist_mat[j][i][idx2] = dist_mat[i][j][idx2]
                if dist_mat[j][i][idx2] != 0:
                    exclude = False
        if exclude:
            continue
        else:
            idx2 += 1

    dist_mat = dist_mat[:idx2]
    if args.mode == 'all':
        with open(args.outputmatrix, 'w') as f:
            f.write(str(idx2) + '\n\n')
            for idx in range(idx2):
                f.write(str(len(tns)) + ' 1 ' + '\n')
                for i in range(len(tns)):
                    f.write(tns[i].label + '     ')
                    for j in range(len(tns)):
                        f.write(str('{:.9f}'.format(dist_mat[i][j][idx])) + ' ')
                    f.write('\n')
                f.write('\n')
    else:
        if args.mode == 'avg':
            summary_dist_mat = np.mean(dist_mat, axis=2)
        elif args.mode == 'min':
            summary_dist_mat = np.min(dist_mat, axis=2)
        elif args.mode == 'med':
            summary_dist_mat = np.median(dist_mat, axis=2)
        else:
            return
        with open(args.outputmatrix, 'w') as f:
            f.write(str(len(tns))+'\n')
            for i in range(len(tns)):
                f.write(tns[i].label + '     ')
                for j in range(len(tns)):
                    f.write(str('{:.9f}'.format(summary_dist_mat[i][j]))+' ')
                f.write('\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CASTLES")
    parser.add_argument("-t", "--speciestree", type=str,  required=True,
                        help="Species tree file in newick format")
    parser.add_argument("-g", "--genetrees", type=str, required=True,
                        help="Gene trees file in newick format")
    parser.add_argument("-m", "--mode", type=str, required=False,
                        help="options: min, avg, med, all", default='avg')
    parser.add_argument("-o", "--outputmatrix", type=str, required=False,
                        help="Output distance matrix in phylip format")
    main(parser.parse_args())
