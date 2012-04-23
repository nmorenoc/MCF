#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm ${dirname}/post_particles.dat 

for F in ${dirname}/mcf_*particles*.out; do
    
    printf "%s" $F
    perl  /scratch/project/trunk/mcf/scripts/particles_3d.perl $F 'post_particles.dat'

done
