#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :                 ppm_module_bem
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains the user-callable subroutines
      !                 of the boundary element solver.
      !                
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_bem.f,v $
      !  Revision 1.4  2004/07/27 09:19:11  oingo
      !  Added ppm_module_bem_quadrule_npoints since it has to be user callable
      !
      !  Revision 1.3  2004/07/26 13:40:28  ivos
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

      MODULE ppm_module_bem

         !----------------------------------------------------------------------
         !  PPM modules
         !----------------------------------------------------------------------
         USE ppm_module_bem_basis
         USE ppm_module_bem_quadrule
         USE ppm_module_bem_quadrule_npoints

      END MODULE ppm_module_bem
