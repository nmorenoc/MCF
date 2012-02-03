#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      ! Module         :            ppm_module_fmm_finalize
      !-------------------------------------------------------------------------
      !
      ! Purpose       : fast multipole method finalize module
      !               
      !
      ! Remarks       :
      !
      ! References    : 
      !
      ! Revisions     :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_fmm_finalize.f,v $
      !  Revision 1.5  2006/06/29 10:28:37  pchatela
      !  Added vector strengths support
      !
      !  Revision 1.3  2005/08/23 14:23:43  polasekb
      !  no difference between single/double
      !
      !  Revision 1.2  2005/05/27 08:42:39  polasekb
      !  removed dummy argument and single/double call
      !
      !  Revision 1.1  2005/05/27 08:03:09  polasekb
      !  initial implementation
      !
      !  Revision 0  2004/11/11 16:35:45  polasekb
      !  Start
      !
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
!#define __SFIELD           9
!#define __VFIELD          10

MODULE ppm_module_fmm_finalize   
  !--------------------------------------------------------------------------
  !Modules
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! INCLUDE THE SOURCES
  !-----------------------------------------------------------------------------

  CONTAINS

#include "ppm_fmm_finalize.F90"

END MODULE ppm_module_fmm_finalize

