#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_fdsolver_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for field solver 
      !                        routines.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
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
#define __SINGLE_PRECISION  1
#define __DOUBLE_PRECISION  2




   MODULE ppm_module_fdsolver_finalize

         !----------------------------------------------------------------------
         !  Define interface to ppm_fdsolver_finalize
         !----------------------------------------------------------------------
         INTERFACE ppm_fdsolver_finalize
            MODULE PROCEDURE ppm_fdsolver_finalize
         END INTERFACE
         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS



#include "ppm_fdsolver_finalize.F90"

      END MODULE ppm_module_fdsolver_finalize
