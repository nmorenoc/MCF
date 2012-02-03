#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_check_topoid
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the utility
      !                 routine which checks topology IDs for their
      !                 validity.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_check_topoid.f,v $
      !  Revision 1.1  2004/08/31 12:13:46  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_check_topoid

         !----------------------------------------------------------------------
         !  Define interface to the topoid check routine
         !----------------------------------------------------------------------
         INTERFACE ppm_check_topoid
            MODULE PROCEDURE ppm_check_topoid
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS
 
#include "ppm_check_topoid.F90"

      END MODULE ppm_module_check_topoid
