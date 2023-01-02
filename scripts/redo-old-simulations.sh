#!/bin/bash

### ASTRAL-II simulations
rep=50
genes=1000
sp=200
for b in 0.000001; do
        for t in 10000000; do
                ./simphy -rs $rep -rl U:$genes,$genes -rg 1 -st U:$t,$t -si U:1,1 -sl U:$sp,$sp -sb U:$b,$b -sp U:200000,200000 -hs ln:1.5,1 -hl ln:1.2,1 -hg ln:1.4,1 -su e:10000000 -so U:1,1 -od 1 -ot 0 -v 3  -cs 293745 -o ../data/ASTRALII/model.$sp.$t.$b 2>&1 1>../data/ASTRALII/model.$sp.$t.$b-simphylog.txt
 		rm  ../data/ASTRALII/model.$sp.$t.$b/$x/g_trees*.trees; 
                rm  ../data/ASTRALII/model.$sp.$t.$b/$x/l_trees.trees
 		rm  ../data/ASTRALII/model.$sp.$t.$b/$x/ultrametric-genetrees.tre;
 		cat ../data/ASTRALII/model.$sp.$t.$b/$x/g_trees*.ralpha > ../data/ASTRALIII/$x/g_trees_ralpha.txt; 
		rm  ../data/ASTRALII/model.$sp.$t.$b/$x/g_trees*.ralpha ; 
        done
done

### MVRoot simulations
g=14.41412

x=(0.15 1.5 5); 
for i in "${x[@]}"; do 
	../ASTRAL-u/scripts/simphy -rs 100 -rl f:500 -rg 1 -sb lu:0.0000001,0.000001 -sd lu:0.0000001,sb -st ln:$g,1 -sl f:30 -si f:1 -sp u:10000,1000000 -su ln:-17.27461,0.6931472 -hh f:1 -hs ln:$i,1 -hl ln:1.551533,0.6931472 -hg ln:$i,1 -cs 9644 -v 3 -o outgroup.0.species.$i.genes.$i -ot 0 -op 1 -od 1 > log.txt; 
	mv log.txt outgroup.0.species.$i.genes.$i/; 
	for x in `seq -w 1 100`; do 
		rm outgroup.0.species.$i.genes.$i/$x/g_trees*.trees; 
		rm outgroup.0.species.$i.genes.$i/$x/l_trees.trees;
		rm outgroup.0.species.$i.genes.$i/$x/ultrametric-genetrees.tre;
		cat outgroup.0.species.$i.genes.$i/$x/g_trees*.ralpha > outgroup.0.species.$i.genes.$i/$x/g_trees_ralpha.txt; 
		rm outgroup.0.species.$i.genes.$i/$x/g_trees*.ralpha ; 
	done
	echo $i; 
done

x=(0.15 1.5 5); 
for i in "${x[@]}"; do 
	../ASTRAL-u/scripts/simphy -rs 100 -rl f:500 -rg 1 -sb lu:0.0000001,0.000001 -sd lu:0.0000001,sb -st ln:$g,1 -sl f:30 -so f:1 -si f:1 -sp u:10000,1000000 -su ln:-17.27461,0.6931472 -hh f:1 -hs ln:$i,1 -hl ln:1.551533,0.6931472 -hg ln:$i,1 -cs 9644 -v 3 -o outgroup.1.species.$i.genes.$i -ot 0 -op 1 -od 1 > log.txt; 
	mv log.txt outgroup.1.species.$i.genes.$i/; 
	echo $i; 
	for x in `seq -w 1 100`; do 
		rm outgroup.1.species.$i.genes.$i/$x/g_trees*.trees; 
		rm outgroup.1.species.$i.genes.$i/$x/l_trees.trees;
		rm outgroup.1.species.$i.genes.$i/$x/ultrametric-genetrees.tre;
		cat outgroup.1.species.$i.genes.$i/$x/g_trees*.ralpha > outgroup.1.species.$i.genes.$i/$x/g_trees_ralpha.txt; 
		rm outgroup.1.species.$i.genes.$i/$x/g_trees*.ralpha ; 
	done
done

exit 
 
x=(0.25 0.5 1 2 4); 
for i in "${x[@]}"; do 
	../ASTRAL-u/scripts/simphy -rs 100 -rl f:500 -rg 1 -sb lu:0.0000001,0.000001 -sd lu:0.0000001,sb -st ln:$g,1 -sl f:30 -so f:$i -si f:1 -sp u:10000,1000000 -su ln:-17.27461,0.6931472 -hh f:1 -hs ln:1.5,1 -hl ln:1.551533,0.6931472 -hg ln:1.5,1 -cs 9644 -v 3 -o outgroup.$i.species.1.5.genes.1.5 -ot 0 -op 1 -od 1 > log.txt; 
	mv log.txt outgroup.$i.species.1.5.genes.1.5/; 
	echo $i;
done

exit 

### ASTRALIII dataset
simphy -rs 50 -rl f:1000 -rg 1 -sb f:0.0000001 -sd f:0 -st ln:14.70055,0.25 -sl f:100 -so f:1 -si f:1 -sp f:400000 -su ln:-17.27461,0.6931472 -hh f:1 -hs ln:1.5,1 -hl ln:1.551533,0.6931472 -hg ln:1.5,1 -cs 9644 -v 3 -o ../data/ASTRALIII -ot 0 -op 1 -od 1 2>&1 1>../data/astral3-simphylog.txt

for x in `seq -w 1 50`; do 
 rm ../data/ASTRALIII/$x/g_trees*.trees; 
 rm ../data/ASTRALIII/$x/l_trees.trees;
 rm ../data/ASTRALIII/$x/ultrametric-genetrees.tre;
 cat ../data/ASTRALIII/$x/g_trees*.ralpha > ../data/ASTRALIII/$x/g_trees_ralpha.txt; 
 rm ../data/ASTRALIII/$x/g_trees*.ralpha ; 
done
