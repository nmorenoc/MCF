#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :           ppm_module_bem_quadrule_npoints
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
      !  $Log: ppm_module_bem_quadrule_npoints.f,v $
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

      MODULE ppm_module_bem_quadrule_npoints

        !-----------------------------------------------------------------------
        !  Define interface to ppm_bem_get_quadrule_points
        !-----------------------------------------------------------------------
        INTERFACE ppm_bem_quadrule_npoints
           MODULE PROCEDURE ppm_bem_quadrule_npoints
        END INTERFACE

        !-----------------------------------------------------------------------
        ! Include the sources
        !-----------------------------------------------------------------------
        CONTAINS

#include "ppm_bem_quadrule_npoints.F90"

      END MODULE ppm_module_bem_quadrule_npoints
