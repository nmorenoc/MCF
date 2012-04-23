#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm ${dirname}/post_particles.dat 

for F in ${dirname}/mcf_*particles*.out; do
    
    printf "%s" $F
    perl  ~/MCF/mcf/scripts/particles_2d.perl $F 'post_particles.dat'

done
