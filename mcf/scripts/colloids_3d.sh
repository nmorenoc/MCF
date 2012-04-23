#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm ${dirname}/colloids.dat 

for F in ${dirname}/mcf_colloid*.out; do
    
    printf "%s" $F
    perl  ~/project/trunk/mcf/scripts/colloids_3d.perl $F ${dirname}/colloids.dat

done
