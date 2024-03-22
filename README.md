# CASTLES

**CASTLES** is a method for estimating branch lengths of a given species tree from estimated gene trees in the unit of expected number of substitutions per sequence site (substitution units), that addresses gene tree heterogeneity due to incomplete lineage sorting (ILS), as modeled by the multi-species coalescent (MSC) model. 

The CASTLES algorithm is described in the following paper:

Y. Tabatabaee, C. Zhang, T. Warnow, S. Mirarab, Phylogenomic branch length estimation using quartets, Bioinformatics, Volume 39, Issue Supplement_1, June 2023, Pages i185–i193, https://doi.org/10.1093/bioinformatics/btad221

Datasets and results from this study are available in [CASTLES-paper](https://github.com/ytabatabaee/CASTLES-paper/tree/main) repository.

**NEW:** An improved version of CASTLES that also handles gene duplication and loss (called **CASTLES-Pro**) is now integerated inside the species tree estimation software [ASTER](https://github.com/chaoszhang/ASTER). We recommend using ASTER directly to get branch lengths on trees produced by [ASTRAL](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral.md), [ASTRAL-Pro](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral-pro.md), or on a fixed input tree topology. 

## Integrated inside ASTER (NEW)
Follow the installation instructions on [ASTER](https://github.com/chaoszhang/ASTER) repository and download ASTER (>= v1.16.2.4). ASTRAL and ASTRAL-Pro with the following compilation by default produce species trees with SU branch lengths. See [ASTRAL tutorial](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral.md) and [ASTRAL-Pro tutorial](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral-pro.md) for more information.
### Compilation
For **single-copy** gene trees, use the following command to compile ASTER
```
$ g++ -D ASTRALIV -std=gnu++11 -march=native -Ofast -pthread src/astral.cpp -o bin/astral4
```
For **multi-copy** gene trees, use the following command for compilation
```
$ g++ -D CASTLES -std=gnu++11 -march=native -Ofast -pthread src/astral-pro.cpp -o bin/astral-pro2
```
### Usage
**Arguments**
```
Required
 -i,  --input        input gene trees in newick format       
 -o,  --output       output species tree with SU branch lengths
Optional
 -i,  --genelength   average gene sequence length [default: 1000]   
 -o,  --root         outgroup name
 -a,  --mapping      list of gene name to taxon name maps
 -c,  --constraint   species tree to score 
```
To infer a species tree using ASTRAL with SU branch lengths, use the following command:
```
$ bin/astral4 -i <gene-tree-path> -o <output-path> [--root <outgroup-name>] [--genelength <gene-length>]
```
To infer branch lengths on a **fixed** species tree topology, use the scoring option `-C -c <species-tree-path>`:
```
$ bin/astral4 -i <gene-tree-path> -C -c <species-tree-path> -o <output-path> [--root <outgroup-name>] [--genelength <gene-length>]
```
To infer branch lengths using ASTRAL-Pro, use the following command:
```
$ bin/astral-pro2 -i <gene-tree-path > [-C -c <species-tree-path>] -o <output-path> [--root <outgroup-name>] [--genelength <gene-length>]
```
If an outgroup is known, it is recommded to specify it using the option `--root`. Additionally, if the average gene sequence length that was used to infer the input gene trees is approximately known, it can be set with `--genelength`. The default value for this parameter is 1000.

#### Handling multiple individuals per species
When there are multiple individuals per species and the individual names do not match the species names, use the following command
```
$ bin/astral4 -i <gene-tree-path> [-C -c <species-tree-path>] -a <name_map> -o <output-path>
```
where the `name_map` file contains maps from individual names to species names in the following format
```
individual_name1    species_name1
individual_name2    species_name2
individual_name3    species_name3
...
```
### Additional Files
- An old documentation of CASTLES is available [here](https://github.com/ytabatabaee/CASTLES/blob/main/OLD-README.md).
- A modified version of the simulation software [SimPhy](https://github.com/adamallo/SimPhy) that produces species trees with SU branch lenghts is available in [simulation_files](https://github.com/ytabatabaee/CASTLES/tree/main/simulation_files).
- Some useful scripts for comparing branch lengths of two species trees are available in [scripts](https://github.com/ytabatabaee/CASTLES/tree/main/scripts).
