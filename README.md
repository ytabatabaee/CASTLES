# CASTLES

**CASTLES** is a method for estimating branch lengths of a given species tree from estimated gene trees in the unit of expected number of substitutions per sequence site (substitution units), that addresses gene tree heterogeneity due to incomplete lineage sorting (ILS), as modeled by the multi-species coalescent (MSC) model.

The CASTLES algorithm is described in the following paper:

Y. Tabatabaee, C. Zhang, T. Warnow, S. Mirarab, Phylogenomic branch length estimation using quartets, Bioinformatics, Volume 39, Issue Supplement_1, June 2023, Pages i185–i193, https://doi.org/10.1093/bioinformatics/btad221

Datasets and results from this study are available in [CASTLES-paper](https://github.com/ytabatabaee/CASTLES-paper/tree/main) repository.

## Dependencies
CASTLES is implemented in Python 3. It was developed and tested in Python version 3.7.0 and has the following dependencies:
- [Python 3.x](https://www.python.org)
- [Dendropy 4.x](https://dendropy.org/index.html)
- [Numpy](https://numpy.org)

## Usage Instructions

**Input:** A file containing a species tree and a file containing a set of single-copy or multi-copy gene trees, both in newick format.

**Output:** A file containing the species tree in newick format, annotated with substitution unit (SU) branch lengths.

Running CASTLES is currently a two-step approach, but in the future it will be integerated inside the species tree estimation software [ASTER](https://github.com/chaoszhang/ASTER) and can be used with both [ASTRAL](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral.md) and [ASTRAL-Pro](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral-pro.md).
1) Annotate branches of the species tree with quartet statistics using [ASTER](https://github.com/chaoszhang/ASTER).
2) Assign final branch lengths to each branch of the species tree using `castles.py`.

### Annotating the input species tree with ASTER
Follow the installation instructions on [ASTER](https://github.com/chaoszhang/ASTER) repository and download ASTER (>= v1.13.2.4). 

For **single-copy** gene trees, use the following command to compile ASTER
```
$ g++ -std=gnu++11 -D"ASTRALIV" -march=native -Ofast -pthread src/astral.cpp -o bin/astral_castles
```
For **multi-copy** gene trees, use the following command for compilation
```
$ g++ -std=gnu++11 -march=native -D CASTLES -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro_castles
```
**WARNING:** *If you use an M1 Mac and you got a compilation error with the command above, we suggest you switch to a Linux system.*

Then, use the following command to run ASTER, where the annotated tree is printed to the log file
```
$ astral_castles -C -i <gene_tree_path> -c <species_tree_path> -o <output_path> > annotated.tre
```
or the following command for **multi-copy** gene trees
```
$ astral-pro_castles -C -i <gene_tree_path> -c <species_tree_path> -o <output_path> > annotated.tre
```
If an outgroup taxon is known, it can be specified with the option `--root` as follows
```
$ astral_castles -C -i <gene_tree_path> -c <species_tree_path> -o <output_path> --root <outgroup_taxon> > annotated.tre
```
#### Handling multiple individuals per species
When there are multiple individuals per species and the individual names do not match the species names, run the following command
```
$ astral_castles -C -i <gene_tree> -a <name_map> -c <species_tree> -o <output_path> > annotated.tre
```
where the `name_map` file contains maps from individual names to species names in the following format
```
individual_name1    species_name1
individual_name2    species_name2
individual_name3    species_name3
...
```
### Assigning SU branch lengths
Use the following command to produce the final species tree with SU branch lengths (**note:** the input is the *ASTER-annotated* tree, not the original species tree)
```
$ python3 castles.py -t annotated.tre -g <gene_tree_path> -o <output_path>
```
**Arguments**
- **Required**
```
 -t,  --speciestree        ASTER-annotated species tree in newick format
 -g,  --genetrees          input single-copy gene trees in newick format
 -o,  --output             output file containing a species tree annotated with SU branch lengths
```


**Example**

The `example` directory contains a 30-taxon model species tree and a corresponding tree annotated by ASTER, and 500 estimated gene trees. The command below shows how CASTLES can be run on this data:
```
$ python3 castles.py -t example/aster.trees.annotated -g example/estimatedgenetre.gtr -o example/castles.tre
```
