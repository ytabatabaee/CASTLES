#!/bin/bash


q=nothing
rm -r $q;
./simphy -rs 200 -rl f:10000 -rg 1 -sb f:0.0001 -sd f:0 -st f:400 -sl f:4 -so f:0 -si f:1 -sp f:100 -su f:0.01 -hs f:10000 -hl f:10000 -hh f:10000 -hg f:10000 -cs 12964 -v 4 -o $q -ot 0 -op 1 -od 1 > log-$q.txt
./postsim.sh $q


q=only_hs
rm -r $q;
./simphy -rs 200 -rl f:10000 -rg 1 -sb f:0.0001 -sd f:0 -st f:400 -sl f:4 -so f:0 -si f:1 -sp f:100 -su ln:-8.27461,0.4 -hs ln:1.5,1 -hl f:10000 -hh f:10000 -hg f:10000 -cs 12964 -v 4 -o $q -ot 0 -op 1 -od 1 > log-$q.txt
./postsim.sh $q


q=only_hl
rm -r $q;
./simphy -rs 200 -rl f:10000 -rg 1 -sb f:0.0001 -sd f:0 -st f:400 -sl f:4 -so f:0 -si f:1 -sp f:100 -su ln:-8.27461,0.4 -hs f:10000 -hl ln:1.551533,0.693 -hh f:10000 -hg f:10000 -cs 12964 -v 4 -o $q -ot 0 -op 1 -od 1 > log-$q.txt
./postsim.sh $q


q=hs_hl
rm -r $q;
./simphy -rs 200 -rl f:10000 -rg 1 -sb f:0.0001 -sd f:0 -st f:400 -sl f:4 -so f:0 -si f:1 -sp f:100 -su ln:-8.27461,0.4 -hs ln:1.5,1 -hl ln:1.551533,0.693 -hh f:10000 -hg f:10000 -cs 12964 -v 4 -o $q -ot 0 -op 1 -od 1 > log-$q.txt
./postsim.sh $q


q=hs_hl_hg
rm -r $q;
./simphy -rs 200 -rl f:10000 -rg 1 -sb f:0.0001 -sd f:0 -st f:400 -sl f:4 -so f:0 -si f:1 -sp f:100 -su ln:-8.27461,0.4 -hs ln:1.5,1 -hl ln:1.551533,0.693 -hh f:10000 -hg ln:1.5,1 -cs 12964 -v 4 -o $q -ot 0 -op 1 -od 1 > log-$q.txt
./postsim.sh $q


q=hs_hl_hg_highr
rm -r $q;
./simphy -rs 200 -rl f:10000 -rg 1 -sb f:0.004 -sd f:0 -st f:400 -sl f:4 -so f:0 -si f:1 -sp f:100 -su ln:-8.27461,0.4 -hs ln:1.5,1 -hl ln:1.551533,0.693 -hh f:10000 -hg ln:1.5,1 -cs 12964 -v 4 -o $q -ot 0 -op 1 -od 1 > log-$q.txt
./postsim.sh $q
