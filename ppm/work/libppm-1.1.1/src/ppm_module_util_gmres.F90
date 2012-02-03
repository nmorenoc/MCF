#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :             ppm_module_util_gmres
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the GMRES
      !                 solver
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_util_gmres.f,v $
      !  Revision 1.1  2006/05/11 10:27:00  pchatela
      !  Initial insertion
      !
      !  
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
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2

      MODULE ppm_module_util_gmres
         INTEGER, PARAMETER   :: ppm_gmres_param_success = 0
         INTEGER, PARAMETER   :: ppm_gmres_param_failure = 1
         INTEGER, PARAMETER   :: ppm_gmres_param_maxiter = 2
               
         INTERFACE ppm_util_gmres_solveupper
            MODULE PROCEDURE ppm_util_gmres_solveupper_s
            MODULE PROCEDURE ppm_util_gmres_solveupper_d
         END INTERFACE
         
         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#define __KIND __SINGLE_PRECISION
#include "ppm_util_gmres_solveupper.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_gmres_solveupper.F90"
#undef __KIND

#define __KIND __SINGLE_PRECISION
#include "ppm_util_gmres.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_util_gmres.F90"
#undef __KIND

      END MODULE ppm_module_util_gmres
