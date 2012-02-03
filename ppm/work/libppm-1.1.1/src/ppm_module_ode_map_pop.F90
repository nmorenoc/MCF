#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_ode_map_pop
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_ode_map_pop
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_ode_map_pop.f,v $
      !  Revision 1.2  2004/07/26 12:00:23  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.1  2004/07/26 07:45:48  michaebe
      !  Procedure modules created in the course of atomization.
      !
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

      MODULE ppm_module_ode_map_pop

        !-----------------------------------------------------------------------
        !  Interfaces
        !-----------------------------------------------------------------------

        INTERFACE ppm_ode_map_pop
           MODULE PROCEDURE ppm_ode_map_pop_s
           MODULE PROCEDURE ppm_ode_map_pop_d
        END INTERFACE
        

      CONTAINS
#define __KIND __SINGLE_PRECISION
#include "ppm_ode_map_pop.F90"
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_ode_map_pop.F90"
#undef  __KIND

      END MODULE ppm_module_ode_map_pop


        
