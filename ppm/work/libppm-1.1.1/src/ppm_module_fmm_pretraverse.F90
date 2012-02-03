#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      ! Module         :            ppm_module_fmm_pretraverse
      !-------------------------------------------------------------------------
      !
      ! Purpose       : module for pretraverse
      !
      ! Remarks       :
      !
      ! References    : 
      !
      ! Revisions     :
      !-------------------------------------------------------------------------
      ! $Log: ppm_module_fmm_pretraverse.f,v $
      ! Revision 1.3  2006/06/29 10:28:37  pchatela
      ! Added vector strengths support
      !
      ! Revision 1.2  2005/09/19 13:03:32  polasekb
      ! code cosmetics
      !
      ! Revision 1.1  2005/05/27 08:04:36  polasekb
      ! initial implementation
      !
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

MODULE ppm_module_fmm_pretraverse   

  !-----------------------------------------------------------------------------
  ! Define Interface
  !-----------------------------------------------------------------------------

  INTERFACE ppm_fmm_pretraverse
	MODULE PROCEDURE ppm_fmm_pretraverse_s
	MODULE PROCEDURE ppm_fmm_pretraverse_d
  END INTERFACE

  !-----------------------------------------------------------------------------
  ! INCLUDE THE SOURCES
  !-----------------------------------------------------------------------------

CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_fmm_pretraverse.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fmm_pretraverse.F90"
#undef __KIND

END MODULE ppm_module_fmm_pretraverse


