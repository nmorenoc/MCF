#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      ! Module         :            ppm_module_fmm_traverse
      !-------------------------------------------------------------------------
      !
      ! Purpose       :  fast multipole method module, tree traversing routine
      !               
      !
      ! Remarks       : 
      !
      ! References    : 
      !
      ! Revisions     :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_fmm_traverse.f,v $
      !  Revision 1.5  2006/06/29 10:28:38  pchatela
      !  Added vector strengths support
      !
      !  Revision 1.4  2005/09/19 13:03:32  polasekb
      !  code cosmetics
      !
      !  Revision 1.3  2005/08/04 16:01:59  polasekb
      !  now really checking whether to use single or double prec.
      !
      !  Revision 1.2  2005/07/27 21:11:07  polasekb
      !  adapted to new subroutine call
      !
      !  Revision 1.1  2005/05/27 07:59:57  polasekb
      !  initial implementation
      !
      !  Revision 0  2004/12/02 15:38:33 polasekb
      !  Start
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
!#define __INTEGER          3
!#define __LOGICAL          4
!#define __2D               7
!#define __3D               8
#define __SFIELD           9
#define __VFIELD          10

MODULE ppm_module_fmm_traverse   

  !-----------------------------------------------------------------------------
  ! Define Interface
  !-----------------------------------------------------------------------------

  INTERFACE ppm_fmm_traverse
	MODULE PROCEDURE ppm_fmm_traverse_s_sf
	MODULE PROCEDURE ppm_fmm_traverse_d_sf
	MODULE PROCEDURE ppm_fmm_traverse_s_vf
	MODULE PROCEDURE ppm_fmm_traverse_d_vf 
  END INTERFACE

  !-----------------------------------------------------------------------------
  ! INCLUDE THE SOURCES
  !-----------------------------------------------------------------------------

CONTAINS

#define __KIND __SINGLE_PRECISION
#define __DIM __SFIELD
#include "ppm_fmm_traverse.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fmm_traverse.F90"
#undef __KIND
#undef __DIM

#define __KIND __SINGLE_PRECISION
#define __DIM __VFIELD
#include "ppm_fmm_traverse.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fmm_traverse.F90"
#undef __KIND
#undef __DIM


END MODULE ppm_module_fmm_traverse

