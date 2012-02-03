#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_fmm
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains the user-callable functions of
      !                 the fmm.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_fmm.f,v $
      !  Revision 1.2  2005/09/19 13:03:31  polasekb
      !  code cosmetics
      !
      !  Revision 1.1  2005/05/27 07:55:59  polasekb
      !  initial implementation
      !
      !  Revision 0  2004/11/16 15:31:00  polasekb
      !  start
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_fmm

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_fmm_init
         USE ppm_module_fmm_finalize
         USE ppm_module_fmm_expansion
         USE ppm_module_fmm_potential
 
      END MODULE ppm_module_fmm
