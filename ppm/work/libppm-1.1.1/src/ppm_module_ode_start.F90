#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_ode_start
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_ode_start
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_ode_start.f,v $
      !  Revision 1.2  2004/07/26 11:21:25  michaebe
      !  added dummy interface
      !
      !  Revision 1.1  2004/07/26 07:45:49  michaebe
      !  Procedure modules created in the course of atomization.
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      MODULE ppm_module_ode_start

        !-----------------------------------------------------
        !  Dummy interface
        !-----------------------------------------------------
        INTERFACE ppm_ode_start
           MODULE PROCEDURE ppm_ode_start
        END INTERFACE

      CONTAINS
#include "ppm_ode_start.F90"

      END MODULE ppm_module_ode_start


        
