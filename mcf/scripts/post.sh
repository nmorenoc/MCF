#! /bin/bash


if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm ${dirname}/particles.dat ${dirname}/kin.dat ${dirname}/mom.dat
for F in ${dirname}/mcf_particles*.out; do
#    ../scripts/convert_output.awk $F >> ${F/.out/.mod} 
    cat ${F/.out/.mod}  >> ${dirname}/particles.dat
    awk '{s+=$4^2+$5^2+$6^2} END{print s/NR}' ${F/.out/.mod} >> ${dirname}/kin.dat
    awk '{x+=$4; y+=$5; z+=$6} END{print x/NR, y/NR, z/NR}' ${F/.out/.mod} >> ${dirname}/mom.dat
    printf "%s\n" $F > "/dev/stderr"
    printf "\n" >> ${dirname}/particles.dat
done
