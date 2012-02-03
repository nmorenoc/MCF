#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef __XLF
@PROCESS NOHOT
#endif
      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_util_fft_backward
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the utility
      !                 routines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_fft_backward.f,v $
      !  Revision 1.2  2004/11/05 18:16:38  michaebe
      !  added xlf compiler directives
      !
      !  Revision 1.1  2004/07/26 07:30:12  ivos
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
#define __SINGLE_PRECISION           1
#define __DOUBLE_PRECISION           2
#define __SINGLE_PRECISION_COMPLEX   5
#define __DOUBLE_PRECISION_COMPLEX   6

      MODULE ppm_module_util_fft_backward

         !----------------------------------------------------------------------
         !  Define interface to ppm_util_fft_backward
         !----------------------------------------------------------------------
         INTERFACE ppm_util_fft_backward
            MODULE PROCEDURE ppm_util_fft_backward_2ds
            MODULE PROCEDURE ppm_util_fft_backward_2dd
            MODULE PROCEDURE ppm_util_fft_backward_2dc
            MODULE PROCEDURE ppm_util_fft_backward_2dcc

            MODULE PROCEDURE ppm_util_fft_backward_3ds
            MODULE PROCEDURE ppm_util_fft_backward_3dd
            MODULE PROCEDURE ppm_util_fft_backward_3dc
            MODULE PROCEDURE ppm_util_fft_backward_3dcc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#define __KIND __SINGLE_PRECISION
#include "ppm_util_fft_backward_2d.F90"
#include "ppm_util_fft_backward_3d.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_fft_backward_2d.F90"
#include "ppm_util_fft_backward_3d.F90"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_util_fft_backward_2d.F90"
#include "ppm_util_fft_backward_3d.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_util_fft_backward_2d.F90"
#include "ppm_util_fft_backward_3d.F90"
#undef __KIND

      END MODULE ppm_module_util_fft_backward
