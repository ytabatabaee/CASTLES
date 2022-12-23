#!/bin/bash

### ASTRALIII dataset
simphy -rs 50 -rl f:1000 -rg 1 -sb f:0.0000001 -sd f:0 -st ln:14.70055,0.25 -sl f:100 -so f:1 -si f:1 -sp f:400000 -su ln:-17.27461,0.6931472 -hh f:1 -hs ln:1.5,1 -hl ln:1.551533,0.6931472 -hg ln:1.5,1 -cs 9644 -v 3 -o ../data/ASTRALIII -ot 0 -op 1 -od 1 2>&1 1>../data/astral3-simphylog.txt

for x in `seq -w 1 50`; do 
 rm ../data/ASTRALIII/$x/g_trees*.trees; 
 rm ../data/ASTRALIII/$x/l_trees.trees;
 rm ../data/ASTRALIII/$x/ultrametric-genetrees.tre;
 cat ../data/ASTRALIII/$x/g_trees*.ralpha > ../data/ASTRALIII/$x/g_trees_ralpha.txt; 
 rm ../data/ASTRALIII/$x/g_trees*.ralpha ; 
done
