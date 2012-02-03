#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ppm_define.h"
#ifdef __XLF
@PROCESS NOOPT
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

#define __INIT             11
#define __SLAB             12


   MODULE ppm_module_fdsolver_solve

         !----------------------------------------------------------------------
         !  Define interface to ppm_fdsolver_solve
         !----------------------------------------------------------------------
 
         INTERFACE ppm_fdsolver_solve_slab
            MODULE PROCEDURE ppm_fdsolver_solve_slab_3d_ss
            MODULE PROCEDURE ppm_fdsolver_solve_slab_3d_sd
            MODULE PROCEDURE ppm_fdsolver_solve_slab_3d_vs
            MODULE PROCEDURE ppm_fdsolver_solve_slab_3d_vd
         END INTERFACE

         INTERFACE ppm_fdsolver_solve_init
            MODULE PROCEDURE ppm_fdsolver_solve_init_2d_ss
            MODULE PROCEDURE ppm_fdsolver_solve_init_2d_sd
            MODULE PROCEDURE ppm_fdsolver_solve_init_2d_vs
            MODULE PROCEDURE ppm_fdsolver_solve_init_2d_vd
            MODULE PROCEDURE ppm_fdsolver_solve_init_3d_ss
            MODULE PROCEDURE ppm_fdsolver_solve_init_3d_sd
            MODULE PROCEDURE ppm_fdsolver_solve_init_3d_vs
            MODULE PROCEDURE ppm_fdsolver_solve_init_3d_vd
         END INTERFACE


         INTERFACE ppm_fdsolver_solve
            MODULE PROCEDURE ppm_fdsolver_solve_2d_ss
            MODULE PROCEDURE ppm_fdsolver_solve_2d_sd
            MODULE PROCEDURE ppm_fdsolver_solve_2d_vs
            MODULE PROCEDURE ppm_fdsolver_solve_2d_vd
            MODULE PROCEDURE ppm_fdsolver_solve_3d_ss
            MODULE PROCEDURE ppm_fdsolver_solve_3d_sd
            MODULE PROCEDURE ppm_fdsolver_solve_3d_vs
            MODULE PROCEDURE ppm_fdsolver_solve_3d_vd
         END INTERFACE
         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __CASE __SLAB

#define __DIM __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND
#undef __DIM



#define __DIM __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND
#undef __DIM

#undef __CASE

#define __CASE __INIT

#define __DIM __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_solve_2d.F90"
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_solve_2d.F90"
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND
#undef __DIM



#define __DIM __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_solve_2d.F90"
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_solve_2d.F90"
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND
#undef __DIM

#undef __CASE




#define __DIM __SFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_solve_2d.F90"
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_solve_2d.F90"
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND
#undef __DIM



#define __DIM __VFIELD
#define __KIND __SINGLE_PRECISION
#include "ppm_fdsolver_solve_2d.F90"
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_fdsolver_solve_2d.F90"
#include "ppm_fdsolver_solve_3d.F90"
#undef  __KIND
#undef __DIM


      END MODULE ppm_module_fdsolver_solve
