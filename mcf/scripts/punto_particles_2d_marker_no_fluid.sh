#! /bin/bash
#############################################################
#This bash script is necessary, as it can be called globally
#after being placed in th ${PATH}.
#The path of perl routine has to be specified mannualy.
#############################################################

#############################################################
#The 1st is radius.
#The 2nd,3rd  parameters are computational box length,
#assuming that x direction is periodic [0,Lx] and
#y direction is wall [0,Ly].
#If a 4th parameter is given, use it as number of colloids,
#otherwise it is zero.
#if a 5th parameter is given, use it as folder name
#otherwise use current directory.
#
#fluid particles are removed.
#############################################################

R=$1
Lx=$2
Ly=$3

if [ -z "$4" ]; then
    Nc=0
else
    Nc=$4
fi
if [ -z "$5" ]; then
    dirname=.
else
    dirname=$5
fi

#############################################################
#remove existing target output file.
#############################################################

rm punto_particles_marker_no_fluid.dat

#############################################################
#start buliding up a punto particles file using perl script.
#############################################################

perl  ~/MCF/mcf/scripts/particles_2d_marker_no_fluid.perl ${R} ${Lx} ${Ly} ${Nc} ${dirname}

