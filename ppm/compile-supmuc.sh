#! /bin/bash

. /etc/profile
. /etc/profile.d/modules.sh

module unload mpi.ibm
module load mpi.intel

make PREFIX=/home/hlrb2/pr32ma/lu64yuz3/MCF/ppm/local_intelmpi_ifort121_o2_g \
     BUILD_MAKE_FLAGS=-j16 FCFLAGS='-g -O2' \
     WRKDIR=/home/hlrb2/pr32ma/lu64yuz3/MCF/ppm/work_intelmpi_ifort121_o2_g \
     install_all_but_mpi
