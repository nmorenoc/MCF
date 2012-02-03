#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      ! Module         :            ppm_module_mg_res
      !-------------------------------------------------------------------------
      !
      ! Purpose       : multigrid module
      !               
      !
      ! Remarks       :
      !
      ! References    : 
      !
      ! Revisions     :
      !-------------------------------------------------------------------------
      ! $Log: ppm_module_mg_res.f,v $
      ! Revision 1.1  2004/09/22 18:45:11  kotsalie
      ! MG new version
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------


#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __INTEGER          3
#define __LOGICAL          4
#define __2D               7
#define __3D               8
#define __SFIELD           9
#define __VFIELD          10

MODULE ppm_module_mg_res  
  !--------------------------------------------------------------------------
  !Modules
  !-----------------------------------------------------------------------------
  
  !---------------------------------------------------------------------------
  INTERFACE ppm_mg_res_sca
     MODULE PROCEDURE ppm_mg_res_coarse_2D_sca_s
     MODULE PROCEDURE ppm_mg_res_coarse_2D_sca_d
     MODULE PROCEDURE ppm_mg_res_coarse_3D_sca_s
     MODULE PROCEDURE ppm_mg_res_coarse_3D_sca_d
     MODULE PROCEDURE ppm_mg_res_fine_2D_sca_s
     MODULE PROCEDURE ppm_mg_res_fine_2D_sca_d
     MODULE PROCEDURE ppm_mg_res_fine_3D_sca_s
     MODULE PROCEDURE ppm_mg_res_fine_3D_sca_d
  END INTERFACE

  INTERFACE ppm_mg_res_vec
     MODULE PROCEDURE ppm_mg_res_coarse_2D_vec_s
     MODULE PROCEDURE ppm_mg_res_coarse_2D_vec_d
     MODULE PROCEDURE ppm_mg_res_coarse_3D_vec_s
     MODULE PROCEDURE ppm_mg_res_coarse_3D_vec_d
     MODULE PROCEDURE ppm_mg_res_fine_2D_vec_s
     MODULE PROCEDURE ppm_mg_res_fine_2D_vec_d
     MODULE PROCEDURE ppm_mg_res_fine_3D_vec_s
     MODULE PROCEDURE ppm_mg_res_fine_3D_vec_d
  END INTERFACE
   


  !-----------------------------------------------------------------------------
  ! INCLUDE THE SOURCES
  !-----------------------------------------------------------------------------

CONTAINS

#define __DIM __SFIELD
#define __MESH_DIM __2D
#define __KIND __SINGLE_PRECISION
#include "ppm_mg_res_coarse.F90"
#include "ppm_mg_res_fine.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_mg_res_coarse.F90"
#include "ppm_mg_res_fine.F90"
#undef __KIND
#undef __MESH_DIM 

#define __MESH_DIM __3D
#define __KIND __SINGLE_PRECISION
#include "ppm_mg_res_coarse.F90"
#include "ppm_mg_res_fine.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_mg_res_coarse.F90"
#include "ppm_mg_res_fine.F90"
#undef __KIND
#undef __MESH_DIM
#undef __DIM

#define __DIM __VFIELD
#define __MESH_DIM __2D
#define __KIND __SINGLE_PRECISION
#include "ppm_mg_res_coarse.F90"
#include "ppm_mg_res_fine.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_mg_res_coarse.F90"
#include "ppm_mg_res_fine.F90"
#undef __KIND
#undef __MESH_DIM

#define __MESH_DIM __3D
#define __KIND __SINGLE_PRECISION
#include "ppm_mg_res_coarse.F90"
#include "ppm_mg_res_fine.F90"
#undef __KIND

#define __KIND __DOUBLE_PRECISION
#include "ppm_mg_res_coarse.F90"
#include "ppm_mg_res_fine.F90"
#undef __KIND
#undef __MESH_DIM
#undef __DIM



END MODULE ppm_module_mg_res


