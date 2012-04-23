#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm ${dirname}/particles_punto_sid.dat 

for F in ${dirname}/mcf_*particles*.out; do
    
    perl  /scratch/project/trunk/mcf/scripts/post_3d_sid.perl $F 'particles_punto_sid.dat'
    printf "%s\n" $F > "/dev/stderr"
    printf "\n" >> ${dirname}/particles_punto_sid.dat

done
