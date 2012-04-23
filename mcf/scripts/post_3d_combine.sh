#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm ${dirname}/particles.dat 

for F in ${dirname}/mcf_particles*.out; do
    
    #printf "%s\n" $F > "/dev/stderr"
    printf "%s%s\n" $F " >> ${dirname}/particles.dat"
    perl  /scratch/project/trunk/mcf/scripts/post_3d.perl $F ${dirname}/particles.dat
done
