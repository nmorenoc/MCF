#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_comp_mesh_fd
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 concerned with computing various operators on
      !                 meshes using finite differences.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_comp_mesh_fd.f,v $
      !  Revision 1.2  2004/07/26 13:36:04  ivos
      !  Commented CONTAINS since NEC and XLF do not like empty modules.
      !
      !  Revision 1.1  2004/07/26 07:29:26  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      !-------------------------------------------------------------------------
      !  Define types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION         1
#define __DOUBLE_PRECISION         2
#define __SINGLE_PRECISION_COMPLEX 3
#define __DOUBLE_PRECISION_COMPLEX 4

      MODULE ppm_module_comp_mesh_fd

         !----------------------------------------------------------------------
         !  Define interfaces to the finite difference routines
         !----------------------------------------------------------------------
!         INTERFACE ppm_comp_mesh_fd
!            MODULE PROCEDURE ppm_comp_mesh_fd_s
!            MODULE PROCEDURE ppm_comp_mesh_fd_d
!            MODULE PROCEDURE ppm_comp_mesh_fd_sc
!            MODULE PROCEDURE ppm_comp_mesh_fd_dc
!         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
!         CONTAINS

#define __KIND __SINGLE_PRECISION
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#undef __KIND

      END MODULE ppm_module_comp_mesh_fd
