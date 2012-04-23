#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm ${dirname}/particles.dat 

for F in ${dirname}/mcf_particles*.out; do
    
    printf "%s" $F
    perl  ~/project/trunk/mcf/scripts/particles_3d.perl $F ${dirname}/particles.dat

done
