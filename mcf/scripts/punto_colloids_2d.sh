#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm ${dirname}/punto_colloids.dat 

for F in ${dirname}/mcf_colloid*.out; do
    
    printf "%s" $F
    perl  ~/MCF/mcf/scripts/colloids_2d.perl $F 'punto_colloids.dat'
    printf "\n" >> ${dirname}/punto_colloids.dat

done
