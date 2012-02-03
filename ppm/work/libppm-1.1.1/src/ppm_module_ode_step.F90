#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :            ppm_module_ode_step
      !-------------------------------------------------------------------------
      !
      !  Purpose      : procedure module for ppm_ode_step
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_ode_step.f,v $
      !  Revision 1.5  2004/07/26 14:58:05  michaebe
      !  renamed the preprocessing defines vector and scalar
      !
      !  Revision 1.4  2004/07/26 12:00:24  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.3  2004/07/26 11:48:28  michaebe
      !  added vector and scalar defines
      !
      !  Revision 1.2  2004/07/26 08:14:15  michaebe
      !  Added overloading for scalar lda.
      !
      !  Revision 1.1  2004/07/26 07:45:50  michaebe
      !  Procedure modules created in the course of atomization.
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Define data types
      !-------------------------------------------------------------------------
#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __SCA 3
#define __VEC 4

      MODULE ppm_module_ode_step

        !-----------------------------------------------------
        !  Interface
        !-----------------------------------------------------

        INTERFACE ppm_ode_step
           MODULE PROCEDURE ppm_ode_step_ss
           MODULE PROCEDURE ppm_ode_step_ds
           MODULE PROCEDURE ppm_ode_step_sv
           MODULE PROCEDURE ppm_ode_step_dv
        END INTERFACE

      CONTAINS
#define __MODE __SCA
#define __KIND __SINGLE_PRECISION
#include "ppm_ode_step.F90"
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_ode_step.F90"
#undef  __KIND
#undef  __MODE
#define __MODE __VEC
#define __KIND __SINGLE_PRECISION
#include "ppm_ode_step.F90"
#undef  __KIND
#define __KIND __DOUBLE_PRECISION
#include "ppm_ode_step.F90"
#undef  __KIND
#undef  __MODE

      END MODULE ppm_module_ode_step


        
