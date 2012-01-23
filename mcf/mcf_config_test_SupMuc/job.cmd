#!/bin/bash
#
#@ job_type = MPICH
#@ class = test
#@ node  = 2
#@ total_tasks=64
#@ wall_clock_limit = 2:00:00
#@ job_name = 3D_64_cores
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = /work/pr32ma/lu64yuz3/scalability/A64/
#@ error = job$(jobid).err
#@ notification=always
#@ notify_user=xin.bian@aer.mw.tum.de
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
module unload mpi.ibm
module load mpi.mpich2/1.4/intel
mpiexec -n 64 ./mcf | tee mcf_std.out
