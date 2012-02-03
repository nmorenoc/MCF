#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module ppm_module_map_part_ghost_test 
contains
  subroutine ppm_module_map_part_ghost_run
#ifdef HAVE_MPI
    use mpi, only : MPI_Init, MPI_COMM_WORLD
#endif
use mod_unit
use ppm_module_map_part_ghost
use ppm_module_mktopo, only : ppm_topo_mkgeom_d
use ppm_module_data, only : ppm_kind_double, ppm_char, ppm_param_decomp_bisection, &
     ppm_param_assign_internal, ppm_param_bcdef_periodic, &
     ppm_param_map_ghost_get, ppm_param_map_push, ppm_param_map_send, ppm_param_map_pop, &
     ppm_param_bcdef_LE
use ppm_module_init, only : ppm_init
use ppm_module_impose_part_bc, only : ppm_impose_part_bc
implicit none

integer, parameter :: MK        = ppm_kind_double
integer, parameter :: MAXCHAR   = ppm_char + 1
real(MK), dimension(:), allocatable                        :: min_phys,max_phys, box_length

integer, parameter  :: ndim    = 2
real(MK)          :: cutoff
integer, parameter :: tolexp = INT(LOG10(EPSILON(cutoff)))+1
integer, dimension(:), allocatable                        :: bcdef 

real(MK), dimension(:  ), pointer  :: sub_cost => null()
integer , dimension(:  ), pointer  :: sub2proc => null()

!> NOTE: fortran 95 features
real(MK), dimension(:,:), pointer  :: min_sub => null() 
real(MK), dimension(:,:), pointer  :: max_sub => null() 
integer , dimension(:  ), pointer  :: isublist => null()
real(MK), dimension(:, :), pointer :: xp => null()

real(MK), dimension(ndim, 2*ndim)  :: shear_length

integer                                       :: nsublist
integer                            :: nsubs
integer                            :: Npart, Mpart

integer, parameter                            :: isymm = 0

integer                            :: topo_id
integer                            :: decomp
integer             :: info

#ifdef HAVE_MPI
call unit_init(5)
call MPI_Init(info)
call unit_assert_equal("MPI_Init", info, 0)
call ppm_init(ndim, MK, tolexp, MPI_COMM_WORLD, 0, info)
call unit_assert_equal("ppm_init", info, 0)
#else 
call unit_init(4)
call ppm_init(ndim, MK, tolexp, 0, 0, info)
call unit_assert_equal("ppm_init", info, 0)
#endif

cutoff = 0.1_MK
decomp = ppm_param_decomp_bisection
allocate(bcdef(2*ndim))
shear_length(1:ndim, 1) = 0.0_MK
shear_length(1:ndim, 2) = 0.0_MK
shear_length(1:ndim, 3) = (/ -0.2_MK, 0.0_MK /)
shear_length(1:ndim, 4) = (/ 0.2_MK, 0.0_MK /)

bcdef(1) = ppm_param_bcdef_periodic
bcdef(2) = ppm_param_bcdef_periodic
bcdef(3) = ppm_param_bcdef_LE
bcdef(4) = ppm_param_bcdef_LE

allocate(max_phys(ndim))
allocate(min_phys(ndim))
allocate(box_length(ndim))
min_phys = (/ 0.0_MK, 0.0_MK /)
max_phys = (/ 1.0_MK, 2.0_MK /)

box_length = max_phys - min_phys;

call ppm_topo_mkgeom_d(decomp, ppm_param_assign_internal, min_phys, max_phys, bcdef,  &
     cutoff,topo_id, min_sub, max_sub, sub_cost, sub2proc, nsubs, isublist, &
     nsublist, info)

Npart = 2
allocate(xp(ndim, Npart))
! low boundary 
xp(1:ndim, 1) = (/0.5_MK, 0.5_MK*cutoff /)

! upper boundary 
xp(1:ndim, 2) = (/0.5_MK, max_phys(2) - 0.5_MK*cutoff /)
call ppm_impose_part_bc(xp, Npart, topo_id, info)
call ppm_map_part_ghost(xp,ndim,Npart,Mpart, isymm ,cutoff, ppm_param_map_ghost_get, info, shear_length)
call ppm_map_part_ghost(xp,ndim,Npart,Mpart,isymm,cutoff,ppm_param_map_send,info)
call ppm_map_part_ghost(xp,ndim,Npart,Mpart,isymm,cutoff,ppm_param_map_pop,info)
call unit_assert_equal("Mpart after mapping", Mpart, 2*Npart)


print *, "last test"
call unit_assert_equal_within("distance between real and ghost", xp(1, 1) - xp(1, Npart+1), shear_length(1, 3), 1e-5_MK)
call unit_assert_equal_within("distance between real and ghost", xp(2, 1) - xp(2, Npart+1), -box_length(2), 1e-5_MK)

call unit_finalize()
end subroutine ppm_module_map_part_ghost_run
end module ppm_module_map_part_ghost_test
