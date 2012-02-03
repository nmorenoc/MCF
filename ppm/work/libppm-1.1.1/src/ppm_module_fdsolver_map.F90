#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef __XLF
@PROCESS NOHOT
#endif
      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_fdsolver
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for field solver 
      !                        routines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
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
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __COMPLEX          3
#define __DOUBLE_COMPLEX   4

#define __SFIELD            9
#define __VFIELD           10



   MODULE ppm_module_fdsolver_map

         !----------------------------------------------------------------------
         !  Define interface to ppm_fdsolver_map
         !----------------------------------------------------------------------
         INTERFACE ppm_fdsolver_map
            MODULE PROCEDURE ppm_fdsolver_map_2d_sca_s
            MODULE PROCEDURE ppm_fdsolver_map_2d_sca_d
            MODULE PROCEDURE ppm_fdsolver_map_2d_sca_c
            MODULE PROCEDURE ppm_fdsolver_map_2d_sca_cc
            MODULE PROCEDURE ppm_fdsolver_map_2d_vec_s
            MODULE PROCEDURE ppm_fdsolver_map_2d_vec_d
            MODULE PROCEDURE ppm_fdsolver_map_2d_vec_c
            MODULE PROCEDURE ppm_fdsolver_map_2d_vec_cc
            MODULE PROCEDURE ppm_fdsolver_map_3d_sca_s
            MODULE PROCEDURE ppm_fdsolver_map_3d_sca_d
            MODULE PROCEDURE ppm_fdsolver_map_3d_sca_c
            MODULE PROCEDURE ppm_fdsolver_map_3d_sca_cc
            MODULE PROCEDURE ppm_fdsolver_map_3d_vec_s
            MODULE PROCEDURE ppm_fdsolver_map_3d_vec_d
            MODULE PROCEDURE ppm_fdsolver_map_3d_vec_c
            MODULE PROCEDURE ppm_fdsolver_map_3d_vec_cc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __DIM __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_map_2d.F90"
#include "ppm_fdsolver_map_3d.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_map_2d.F90"
#include "ppm_fdsolver_map_3d.F90"
#undef  __KIND

#define __KIND __COMPLEX
#include "ppm_fdsolver_map_2d.F90"
#include "ppm_fdsolver_map_3d.F90"
#undef  __KIND

#define __KIND __DOUBLE_COMPLEX
#include "ppm_fdsolver_map_2d.F90"
#include "ppm_fdsolver_map_3d.F90"
#undef  __KIND
#undef  __DIM


#define __DIM __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_map_2d.F90"
#include "ppm_fdsolver_map_3d.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_map_2d.F90"
#include "ppm_fdsolver_map_3d.F90"
#undef  __KIND

#define __KIND __COMPLEX
#include "ppm_fdsolver_map_2d.F90"
#include "ppm_fdsolver_map_3d.F90"
#undef  __KIND

#define __KIND __DOUBLE_COMPLEX
#include "ppm_fdsolver_map_2d.F90"
#include "ppm_fdsolver_map_3d.F90"
#undef  __KIND
#undef  __DIM

      END MODULE ppm_module_fdsolver_map
