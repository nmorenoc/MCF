#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm punto_particles_marker_no_fluid.dat

perl  ~/project/trunk/mcf/scripts/particles_3d_marker_no_fluid.perl ${dirname}

