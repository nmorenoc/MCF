#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_decomp_cartesian
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
      !  $Log: ppm_module_decomp_cartesian.f,v $
      !  Revision 1.1  2004/07/26 07:29:30  ivos
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
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_decomp_cartesian

         !----------------------------------------------------------------------
         !  Define interface to cartesian decomposition routine
         !----------------------------------------------------------------------
         INTERFACE ppm_decomp_cartesian
            MODULE PROCEDURE ppm_decomp_cartesian_s
            MODULE PROCEDURE ppm_decomp_cartesian_d
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_decomp_cartesian.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_decomp_cartesian.F90"
#undef __KIND

      END MODULE ppm_module_decomp_cartesian
