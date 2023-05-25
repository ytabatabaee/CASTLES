
### Compare branch lengths of two trees
The script `compare_trees_bl.py` compares the branch lengths of two trees, and outputs several error metrics such as RMSE, mean logarithmic error, etc. 

**Input:** Two phylogenetic trees (with potentially different topologies) in newick format.

**Output:** Several error metrics, branch lengths tabulated in a csv file, and a correlation plot (optional, with `-p`).

```
$ python3 compare_trees_bl.py -t1 <tree1-path> -t2 <tree2-path> [-p]
```
**Example:**
```
$ python3 compare_trees_bl.py -t1 ../example/castles.tre -t2 ../example/s_tree.trees -p
```
**Output:**
```
   Taxon Branch Type        l1        l2
0     23    terminal  0.006886  0.000709
1     27    terminal  0.008258  0.000363
                    ...
56    21    terminal  0.133469  0.138942
57          internal  0.000000  0.000000

Bias: 0.00012074072409021145
Mean absolute error: 0.007556294339972311
Root mean square error (RMSE): 0.010395148343365083
Mean logarithmic error: 0.2767746431429533
Number of negative branches in t2: 0
```

### Compute a patristic distance matrix
The script `patristic_dist_matrix.py` computes a patristic (path-length) distance matrix from a given set of gene trees. The option `-m` specifies the type of distance matrix: `avg`,
`min` and `med` compute a single distance matrix corresponding to the average, minimum and median path-length
distances between pairs of nodes in a set of gene trees. Option `-m all` computes one patristic distance matrix per gene. *Note:* gene trees with all branches zero are not included in distance matrix calculation.

**Input:** A file containing a set of gene trees in newick format and the type of distance matrix specified with `-m`.

**Output:** A file containing a distance matrix in Phylip format.

```
$ python3 patristic_dist_matrix.py -g <gene_tree_path> -m <type> -o <output_path>
```
**Example:**
```
$ python3 patristic_dist_matrix.py -g ../example/estimatedgenetre.gtr -m avg -o ../example/avg_dist_mat.phylip
```
