#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module m_vizing
  implicit none
  private
  public vizing_coloring_wr
  interface
     subroutine vizing_coloring_local(ppm_nproc, &
          nlinks, ilinks, optres) bind(c, name="vizing_coloring_")
       use iso_c_binding, only : c_ptr, c_int
       implicit none
       integer(c_int), intent(inout) :: ppm_nproc
       integer(c_int), intent(inout) :: nlinks
       integer(c_int), dimension(*), intent(inout)  :: ilinks
       integer(c_int), dimension(*), intent(inout)  :: optres
     end subroutine vizing_coloring_local
  end interface
contains
  subroutine vizing_coloring_wr(ppm_nproc, &
       nlinks, ilinks, optres) 
    implicit none
    integer, intent(inout) :: ppm_nproc
    integer, intent(inout) :: nlinks
    integer, dimension(:)  :: ilinks
    integer, dimension(:)  :: optres
    !print *, __FILE__, ':', __LINE__, 'calling vizing_coloring'
    call vizing_coloring_local(ppm_nproc, nlinks, ilinks, optres)
  end subroutine vizing_coloring_wr
end module m_vizing
