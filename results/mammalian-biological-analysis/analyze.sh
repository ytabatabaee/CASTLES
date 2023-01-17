#!/bin/bash

nw_prune ../logs/genetrees.tre GAL > genetees-nogal.tre
python3 ~/Research/projects/local/ASTRAL-u/ASTRAL-u/castles.py -t <(  ~/workspace/ASTER/bin/astral -C -i genetees-nogal.tre -c <( nw_reroot astralv5.7.8.tre GAL ) -o aster-noGALgenetrees.tre --root GAL|nw_prune - GAL  ) -g genetees-nogal.tre  -o castle-test4.tre
