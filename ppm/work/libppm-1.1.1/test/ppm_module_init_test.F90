#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module ppm_module_init_test
contains
  subroutine ppm_module_init_run
    use mod_unit
#ifdef HAVE_MPI
    use mpi, only : MPI_Init, MPI_COMM_WORLD
#endif
     use ppm_module_data, only : ppm_kind_double, ppm_char
     use ppm_module_init, only : ppm_init
     use ppm_module_finalize, only : ppm_finalize
    implicit none
    integer, parameter  :: MK        = ppm_kind_double
    integer, parameter  :: MAXCHAR   = ppm_char + 1
    
    integer, parameter  :: dim    = 2
    real(MK), parameter :: cutoff = 1.0_MK
    integer, parameter  :: tolexp = INT(LOG10(EPSILON(cutoff)))+1
    integer             :: info
    
    info = 0 
#ifdef HAVE_MPI

    call unit_init(4)
    call MPI_Init(info)
    print *, 'preved'
    call unit_assert_equal("MPI_Init", info, 0)


    call ppm_init(dim, MK, tolexp, MPI_COMM_WORLD, 0, info)
    call unit_assert_equal("ppm_init", info, 0)

    call ppm_finalize(info)
    call unit_assert_equal("ppm_finalize", info, 0)

    call MPI_Finalize(info)
    call unit_assert_equal("MPI_Finalize", info, 0)
#else 
    call unit_init(1)
    call ppm_init(dim, MK, tolexp, 0, 0, info)
    call unit_assert_equal("MPI_Init", info, 0)
#endif
    call unit_finalize()
  end subroutine ppm_module_init_run
end module ppm_module_init_test
