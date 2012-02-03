#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_define_subs_bc
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
      !  $Log: ppm_module_define_subs_bc.f,v $
      !  Revision 1.1  2004/07/26 08:55:32  ivos
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

      MODULE ppm_module_define_subs_bc

         !----------------------------------------------------------------------
         !  Define interface to boundary condition definition routine
         !----------------------------------------------------------------------
         INTERFACE ppm_define_subs_bc
            MODULE PROCEDURE ppm_define_subs_bc_s
            MODULE PROCEDURE ppm_define_subs_bc_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_define_subs_bc.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_define_subs_bc.F90"
#undef __KIND

      END MODULE ppm_module_define_subs_bc
