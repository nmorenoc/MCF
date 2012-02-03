#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_mg
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains the user-callable functions of
      !                 the mg solver.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_mg.f,v $
      !  Revision 1.1  2004/09/22 18:39:26  kotsalie
      !  MG new version
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_mg

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_mg_init
         USE ppm_module_mg_solv
         USE ppm_module_mg_finalize
         
      END MODULE ppm_module_mg
