#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :               ppm_module_fieldsolver
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains the user-callable functions of
      !                 the FFT-based fieldsolver.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_fieldsolver.f,v $
      !  Revision 1.3  2004/07/26 13:40:29  ivos
      !  Initial implementation. These are meta-modules for the user-
      !  callable functions. Only these mod files will be given away
      !  to the user.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_fieldsolver

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_fdsolver_poisson
         USE ppm_module_fdsolver_solve
         
      END MODULE ppm_module_fieldsolver
