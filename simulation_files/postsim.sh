#!/bin/bash

d=$1
echo working on $d .....

for x in $d/???; do cat $x/g_trees*trees|sed -e "s/_0_0//g" >> $x/truegenetrees; rm $x/g_trees*trees;  done; 
echo ASTRAL ...;
for x in $d/???; do  java -jar /Users/smirarab/workspace/ASTRALmaster/astral.5.17.4.jar -i $x/truegenetrees -q $x/s_tree.trees -u -o $x/s_tree.trees-estimatedbl-minus2; cat $x/s_tree.trees; done 1>astralrun-$d.txt 2>&1; 
for x in $d/???; do ./get-bl.sh $x; done >  bl-stats-$d.txt; 
grep "7-node species tree correctly simulated" log-$d.txt |sed -e "s/.* //g"  >  $d/species_trees.cu
grep "Global subs" log-$d.txt |tail -n+2 >  $d/mu.txt
grep "" $d/*/s_tree.ralpha |sed -e "s/:/\t/" > $d/branchmu.txt
