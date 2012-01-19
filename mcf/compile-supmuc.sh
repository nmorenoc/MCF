#! /bin/bash

set +e
set +u
source /etc/profile
source /etc/profile.d/modules.sh

module unload mpi.ibm
module load mpi.intel

PREFIX=${HOME}/MCF/ppm/local_intelmpi_ifort121_o2_g

./configure --prefix=$HOME/MCF/mcf/mcf_install/ FC=ifort MPIFC=mpif90 \
    LDFLAGS=-L${PREFIX}/lib/ FCFLAGS="-I${PREFIX}/include -g" \
    MAKEDEPF90=${PREFIX}/bin/makedepf90 

make -j 8