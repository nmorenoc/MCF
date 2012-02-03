#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef __XLF
@PROCESS NOHOT
#endif
      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_util_fft_forward
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
      !  $Log: ppm_module_util_fft_forward.f,v $
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

      MODULE ppm_module_util_fft_forward

         !----------------------------------------------------------------------
         !  Define interface to ppm_util_fft_forward
         !----------------------------------------------------------------------
         INTERFACE ppm_util_fft_forward
            MODULE PROCEDURE ppm_util_fft_forward_2ds
            MODULE PROCEDURE ppm_util_fft_forward_2dd
            MODULE PROCEDURE ppm_util_fft_forward_2dc
            MODULE PROCEDURE ppm_util_fft_forward_2dcc

            MODULE PROCEDURE ppm_util_fft_forward_3ds
            MODULE PROCEDURE ppm_util_fft_forward_3dd
            MODULE PROCEDURE ppm_util_fft_forward_3dc
            MODULE PROCEDURE ppm_util_fft_forward_3dcc
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#define __KIND __SINGLE_PRECISION
#include "ppm_util_fft_forward_2d.F90"
#include "ppm_util_fft_forward_3d.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_fft_forward_2d.F90"
#include "ppm_util_fft_forward_3d.F90"
#undef __KIND

#define __KIND __SINGLE_PRECISION_COMPLEX
#include "ppm_util_fft_forward_2d.F90"
#include "ppm_util_fft_forward_3d.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION_COMPLEX
#include "ppm_util_fft_forward_2d.F90"
#include "ppm_util_fft_forward_3d.F90"
#undef __KIND

      END MODULE ppm_module_util_fft_forward
