#!/bin/bash

d=$1 

echo CU branch lengths:
cat $d/species_trees.cu
echo rates:
paste <( cat $d/mu.txt ) <( cat $d/branchmu.txt )
echo SU branch lengths:
grep "" $d/*/s_tree.trees| sed -e "s/:/\t/"
echo Internal branch lengths:
cat $d/*/s_tree.trees|nw_reroot - 1|nw_reroot -d -|nw_distance -si -mp -;
echo summary:
cat bl-stats-$d.txt |awk '{p=($3-$4/2)/($3+$4);d=-log(1-p);print("p=",p,"d=",d,"e=",($7-$8)*d/(d-p)*(1+2*p)/3,"e2=",($7*(1+2*p)-$8*(1-p))/3,"-----",$0)}'|grep inter|column -t
