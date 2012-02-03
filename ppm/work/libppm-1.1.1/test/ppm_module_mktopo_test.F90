module ppm_module_mktopo_test

contains 

  subroutine ppm_module_mktopo_run
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_MPI
    use mpi, only : MPI_Init, MPI_COMM_WORLD
#endif
    use mod_unit
    use ppm_module_mktopo
    use ppm_module_data, only : ppm_kind_double, ppm_char, ppm_param_decomp_bisection, &
         ppm_param_assign_internal, ppm_param_bcdef_periodic
    use ppm_module_init, only : ppm_init
    implicit none

    integer, parameter :: MK        = ppm_kind_double
    integer, parameter :: MAXCHAR   = ppm_char + 1
    real(MK), dimension(:), allocatable                        :: min_phys,max_phys

    integer, parameter  :: ndim    = 2
    real(MK), parameter :: cutoff = 1.0_MK
    integer, parameter :: tolexp = INT(LOG10(EPSILON(cutoff)))+1
    integer, dimension(:), allocatable                        :: bcdef 

    real(MK), dimension(:  ), pointer  :: sub_cost => null()
    integer , dimension(:  ), pointer  :: sub2proc => null()

    !> NOTE: fortran 95 features
    real(MK), dimension(:,:), pointer  :: min_sub => null() 
    real(MK), dimension(:,:), pointer  :: max_sub => null() 
    integer , dimension(:  ), pointer  :: isublist => null()
    !>  Number of subdomains on local processor
    integer                                       :: nsublist
    integer                            :: nsubs


    integer                            :: topo_id
    integer                            :: decomp
    integer             :: info

    
#ifdef HAVE_MPI
    call unit_init(2)
    call MPI_Init(info)
    call unit_assert_equal("MPI_Init", info, 0)
    call ppm_init(ndim, MK, tolexp, MPI_COMM_WORLD, 0, info)
    call unit_assert_equal("ppm_init", info, 0)
#else 
    call unit_init(1)
    call ppm_init(ndim, MK, tolexp, 0, 0, info)
    call unit_assert_equal("ppm_init", info, 0)
#endif

    decomp = ppm_param_decomp_bisection
    allocate(bcdef(2*ndim))
    bcdef = ppm_param_bcdef_periodic
    ! set domain box
    allocate(max_phys(ndim))
    allocate(min_phys(ndim))
    min_phys = (/ 0.0_MK, 0.0_MK /)
    max_phys = (/ 1.0_MK, 2.0_MK /)

    call ppm_topo_mkgeom_d(decomp, ppm_param_assign_internal, min_phys, max_phys, bcdef,  &
         cutoff,topo_id, min_sub, max_sub, sub_cost, sub2proc, nsubs, isublist, &
         nsublist, info)
    call unit_finalize()
  end subroutine ppm_module_mktopo_run
end module ppm_module_mktopo_test
