# CASTLES

**CASTLES** is a quartet-based method for phylogenomic branch length estimation (in substitution units), that is designed based on the properties of the multispecies coalescent (MSC) model. 

## Dependencies
CASTLES is implemented in Python 3. It was developed and tested in Python version 3.7.0 and has the following dependencies:
- [Python 3.x](https://www.python.org)
- [Dendropy 4.x](https://dendropy.org/index.html)
- [Numpy](https://numpy.org)

## Usage Instructions

**Input:** A file containing a species tree and a file containing a set of gene trees, both in newick format (with or without branch lengths).

**Output:** A file containing the species tree in newick format, annotated with substituion unit (SU) branch lengths.
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
