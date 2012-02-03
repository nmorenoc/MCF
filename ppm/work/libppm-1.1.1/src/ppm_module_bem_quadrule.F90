#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_bem_quadrule
      !-------------------------------------------------------------------------
      !
      !  Purpose      : bem module
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_bem_quadrule.f,v $
      !  Revision 1.1  2004/07/26 07:29:25  ivos
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
      !  Define data types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_bem_quadrule

        !-----------------------------------------------------------------------
        !  Define interface to ppm_bem_quadrule
        !-----------------------------------------------------------------------
        INTERFACE ppm_bem_quadrule
           MODULE PROCEDURE ppm_bem_quadrule_s
           MODULE PROCEDURE ppm_bem_quadrule_d
        END INTERFACE

        !-----------------------------------------------------------------------
        ! Include the sources
        !-----------------------------------------------------------------------
        CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_bem_quadrule.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_bem_quadrule.F90"
#undef __KIND

      END MODULE ppm_module_bem_quadrule
