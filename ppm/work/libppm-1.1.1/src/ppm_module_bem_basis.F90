#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_bem_basis
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
      !  $Log: ppm_module_bem_basis.f,v $
      !  Revision 1.1  2004/07/26 07:29:24  ivos
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

      MODULE ppm_module_bem_basis

        !-----------------------------------------------------------------------
        !  Define interface to ppm_bem_basis
        !-----------------------------------------------------------------------
        INTERFACE ppm_bem_basis
           MODULE PROCEDURE ppm_bem_basis_s
           MODULE PROCEDURE ppm_bem_basis_d
        END INTERFACE

        !-----------------------------------------------------------------------
        ! Include the sources
        !-----------------------------------------------------------------------
        CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_bem_basis.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_bem_basis.F90"
#undef __KIND

      END MODULE ppm_module_bem_basis
