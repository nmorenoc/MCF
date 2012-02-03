#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_ode_alldone
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_ode_allone
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_ode_alldone.f,v $
      !  Revision 1.4  2004/07/26 11:57:45  michaebe
      !  replace module function by module procedure ;)
      !
      !  Revision 1.3  2004/07/26 11:39:49  michaebe
      !  replaced module procedure by module function
      !
      !  Revision 1.2  2004/07/26 11:20:13  michaebe
      !  added dummy interface
      !
      !  Revision 1.1  2004/07/26 07:45:46  michaebe
      !  Procedure modules created in the course of atomization.
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      MODULE ppm_module_ode_alldone

        
      CONTAINS
#include "ppm_ode_alldone.F90"

      END MODULE ppm_module_ode_alldone


        
