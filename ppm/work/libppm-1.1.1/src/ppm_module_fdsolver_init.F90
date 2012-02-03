#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef __XLF
@PROCESS NOHOT
#endif
      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_fdsolver_solve
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
#define __SINGLE_PRECISION  1
#define __DOUBLE_PRECISION  2
#define __COMPLEX           3
#define __DOUBLE_COMPLEX    4

#define __SFIELD            9
#define __VFIELD           10

#define __2D               22
#define __3D               33


   MODULE ppm_module_fdsolver_init

         !----------------------------------------------------------------------
         !  Define interface to ppm_fdsolver_init
         !----------------------------------------------------------------------
         INTERFACE ppm_fdsolver_init
            MODULE PROCEDURE ppm_fdsolver_init_2d_sca_s
            MODULE PROCEDURE ppm_fdsolver_init_2d_sca_d
            MODULE PROCEDURE ppm_fdsolver_init_2d_vec_s
            MODULE PROCEDURE ppm_fdsolver_init_2d_vec_d
            MODULE PROCEDURE ppm_fdsolver_init_3d_sca_s
            MODULE PROCEDURE ppm_fdsolver_init_3d_sca_d
            MODULE PROCEDURE ppm_fdsolver_init_3d_vec_s
            MODULE PROCEDURE ppm_fdsolver_init_3d_vec_d
         END INTERFACE
         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __MESH_DIM __2D
#define __DIM __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_init.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_init.F90"
#undef  __KIND
#undef __DIM


#define __DIM __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_init.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_init.F90"
#undef  __KIND
#undef __DIM
#undef __MESH_DIM




#define __MESH_DIM __3D
#define __DIM __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_init.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_init.F90"
#undef  __KIND
#undef __DIM


#define __DIM __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_init.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_init.F90"
#undef  __KIND
#undef __DIM
#undef __MESH_DIM


      END MODULE ppm_module_fdsolver_init
