#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_ode_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_ode_init
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_ode_init.f,v $
      !  Revision 1.2  2004/07/26 11:21:24  michaebe
      !  added dummy interface
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
      MODULE ppm_module_ode_init

        !-----------------------------------------------------
        !  Dummy interface
        !-----------------------------------------------------
        INTERFACE ppm_ode_init
           MODULE PROCEDURE ppm_ode_init
        END INTERFACE

      CONTAINS
#include "ppm_ode_init.F90"

      END MODULE ppm_module_ode_init


        
