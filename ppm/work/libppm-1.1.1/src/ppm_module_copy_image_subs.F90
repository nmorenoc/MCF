#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_copy_image_subs
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the
      !                 decomposition routines. 
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_copy_image_subs.f,v $
      !  Revision 1.1  2004/07/26 08:55:31  ivos
      !  Renamed.
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

      MODULE ppm_module_copy_image_subs

         !----------------------------------------------------------------------
         !  Define interface to periodic image copy routine
         !----------------------------------------------------------------------
         INTERFACE ppm_copy_image_subs
            MODULE PROCEDURE ppm_copy_image_subs_s
            MODULE PROCEDURE ppm_copy_image_subs_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_copy_image_subs.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_copy_image_subs.F90"
#undef __KIND

      END MODULE ppm_module_copy_image_subs
