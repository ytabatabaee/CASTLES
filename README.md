# CASTLES

**CASTLES** is a quartet-based method for phylogenomic branch length estimation (in substitution units), that is designed based on the properties of the multispecies coalescent (MSC) model.

## Dependencies
CASTLES is implemented in Python 3. It was developed and tested in Python version 3.7.0 and has the following dependencies:
- [Python 3.x](https://www.python.org)
- [Dendropy 4.x](https://dendropy.org/index.html)
- [Numpy](https://numpy.org)

CASTLES must be run on a tree annotated by [ASTER](https://github.com/chaoszhang/ASTER), and in future it will be available to use inside ASTER.

## Usage Instructions

**Input:** A file containing an ASTER-annotated species tree and a file containing a set of gene trees, both in newick format.

**Output:** A file containing the species tree in newick format, annotated with substitution unit (SU) branch lengths.
```
$ python3 castles.py -t <annotated_species_tree_path> -g <gene_tree_path> -o <output_path>
```
**Arguments**
- **Required**
```
 -t,  --speciestree        input species tree in newick format
 -g,  --genetrees          input gene trees in newick format
 -o,  --output             output file containing a species tree annotated with SU branch lengths
```

**Example**
The `example` directory contains one example set, with a 30-taxon model species tree and a corresponding tree annotated by ASTER, and 500 estimated gene trees. The command below shows how CASTLES can be run on this data:
```
$ python3 castles.py -t example/aster.trees.annotated -g example/estimatedgenetre.gtr -o example/castles.tre
```

## Data Availability
Datasets and results from this paper are available in [CASTLES-paper](https://github.com/ytabatabaee/CASTLES-paper/tree/main) repository.
