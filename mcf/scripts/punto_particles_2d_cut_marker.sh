#! /bin/bash

if [ -z "$1" ]; then
    dirname=.
else
    dirname=$1
fi

rm punto_particles_2d_cut_marker.dat

perl  ~/MCF/mcf/scripts/particles_2d_cut_marker.perl ${dirname}

