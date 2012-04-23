#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm punto_particles_pid_marker.dat

perl  ~/MCF/mcf/scripts/particles_pid_2d_marker.perl ${dirname}

