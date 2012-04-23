#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm ${dirname}/particles_punto.dat 

for F in ${dirname}/mcf_*particles*.out; do
    
    #printf "%s\n" $F > "/dev/stderr"
    printf "%s%s\n" $F " >> ${dirname}/particles_punto.dat"
    perl  /scratch/project/trunk/mcf/scripts/post_3d.perl $F ${dirname}/particles_punto.dat
    printf "\n" >> ${dirname}/particles_punto.dat

done
