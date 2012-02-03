#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

       !------------------------------------------------------------------------
       ! Module         :            ppm_module_data_mg
       !------------------------------------------------------------------------
       !
       ! Purpose       : multigrid data module
       !               
       !
       ! Remarks       :
       !
       ! References    : 
       !
       ! Revisions     :
       !------------------------------------------------------------------------
       !  $Log: ppm_module_data_mg.f,v $
       !  Revision 1.6  2006/07/21 11:30:57  kotsalie
       !  FRIDAY
       !
       !  Revision 1.5  2005/12/08 12:42:36  kotsalie
       !  commiting dirichlet
       !
       !  Revision 1.4  2004/10/29 16:00:47  kotsalie
       !  RED BLACK SOR
       !
       !  Revision 1.3  2004/09/28 14:18:19  kotsalie
       !  Added 4th order
       !
       !  Revision 1.2  2004/09/22 18:40:26  kotsalie
       !  MG new version
       !
       !------------------------------------------------------------------------
       !  Parallel Particle Mesh Library (PPM)
       !  Institute of Computational Science
       !  ETH Zentrum, Hirschengraben 84
       !  CH-8092 Zurich, Switzerland
       !------------------------------------------------------------------------


#define __SINGLE_PRECISION 1
#define __DOUBLE_PRECISION 2
#define __INTEGER          3
#define __LOGICAL          4
#define __2D               7
#define __3D               8
#define __SFIELD           9
#define __VFIELD          10 

MODULE ppm_module_data_mg   
   !----------------------------------------------------------------------------
   !Modules
   !----------------------------------------------------------------------------
    USE ppm_module_data,ONLY:ppm_kind_single,ppm_kind_double
    PRIVATE :: ppm_kind_single,ppm_kind_double
   !----------------------------------------------------------------------------
   !The boundary condition values!!
   !----------------------------------------------------------------------------

#define __DIM __SFIELD
#define __MESH_DIM __2D
#define __KIND __SINGLE_PRECISION
   TYPE bc_value_2d_sca_s
      ! 1st index mesh position locally
      REAL(ppm_kind_single), DIMENSION(:),POINTER   ::  pbcvalue
   END TYPE bc_value_2d_sca_s
#undef  __KIND
#define __KIND == __DOUBLE_PRECISION
   TYPE bc_value_2d_sca_d
      !1st index mesh position locally
      REAL(ppm_kind_single), DIMENSION(:),POINTER   ::   pbcvalue
   END TYPE bc_value_2d_sca_d
#undef __KIND



#define __KIND == __SINGLE_PRECISION
   !----------------------------------------------------------------------------
   ! Our multigrid field with all its necessary components (Take a look at the 
   ! theory)
   !----------------------------------------------------------------------------
   TYPE mg_field_2d_sca_s
      !function corrections, error restrictions, errors
      !1st and 2nd index: mesh position(local)
      REAL(ppm_kind_single), DIMENSION(:,:),POINTER  ::  uc 
      REAL(ppm_kind_single), DIMENSION(:,:),POINTER  ::  fc 
      REAL(ppm_kind_single), DIMENSION(:,:),POINTER  ::  err  

      !lets save the boundary condition.index:face of the subdomain(1:4)
      TYPE(bc_value_2d_sca_s), DIMENSION(:), POINTER   ::  bcvalue   
   END TYPE mg_field_2d_sca_s
#undef  __KIND

#define __KIND  == __DOUBLE_PRECISION  

   TYPE mg_field_2d_sca_d
      !function corrections, error restrictions, errors,
      !1st:3rd index: mesh position(local)
      REAL(ppm_kind_double), DIMENSION(:,:),POINTER  ::  uc  
      REAL(ppm_kind_double), DIMENSION(:,:),POINTER  ::  fc  
      REAL(ppm_kind_double), DIMENSION(:,:),POINTER  ::  err  
      !lets save the boundary condition.index:face of the subdomain(1:4)
      TYPE(bc_value_2d_sca_d), DIMENSION(:), POINTER   ::  bcvalue   
   END TYPE mg_field_2d_sca_d
#undef  __KIND


#define __KIND == __SINGLE_PRECISION
   !1st index: subdomain,2nd index : multigrid level
   TYPE(mg_field_2d_sca_s), DIMENSION(:,:),              POINTER  ::   mgfield_2d_sca_s 
#undef __KIND
#define __KIND == __DOUBLE_PRECISION
   !1st index: subdomain,2nd index : multigrid level
   TYPE(mg_field_2d_sca_d), DIMENSION(:,:),              POINTER  ::   mgfield_2d_sca_d 
#undef __KIND
#undef __MESH_DIM

#define __MESH_DIM __3D
#define __KIND __SINGLE_PRECISION
   TYPE bc_value_3d_sca_s
      ! 1st-2nd index mesh position locally
      REAL(ppm_kind_single), DIMENSION(:,:),POINTER   ::  pbcvalue
   END TYPE bc_value_3d_sca_s
#undef  __KIND
#define __KIND == __DOUBLE_PRECISION
   TYPE bc_value_3d_sca_d
      !1st-2nd index mesh position locally
      REAL(ppm_kind_single), DIMENSION(:,:),POINTER   ::   pbcvalue
   END TYPE bc_value_3d_sca_d
#undef __KIND


#define __KIND == __SINGLE_PRECISION
   !----------------------------------------------------------------------------
   ! Our multigrid field with all its necessary components (Take a look at the 
   ! theory)
   !----------------------------------------------------------------------------
   TYPE mg_field_3d_sca_s
      !function corrections, error restrictions, errors
      !1st-3rd index: mesh position(local)
      REAL(ppm_kind_single), DIMENSION(:,:,:),POINTER  ::  uc 
      REAL(ppm_kind_single), DIMENSION(:,:,:),POINTER  ::  fc 
      REAL(ppm_kind_single), DIMENSION(:,:,:),POINTER  ::  err  

      !lets save the boundary condition.index:face of the subdomain(1:6)
      TYPE(bc_value_3d_sca_s), DIMENSION(:), POINTER   ::  bcvalue   
   END TYPE mg_field_3d_sca_s
#undef  __KIND

#define __KIND  == __DOUBLE_PRECISION  

   TYPE mg_field_3d_sca_d
      !function corrections, error restrictions, errors,
      !1st:3rd index: mesh position(local)
      REAL(ppm_kind_double), DIMENSION(:,:,:),POINTER  ::  uc  
      REAL(ppm_kind_double), DIMENSION(:,:,:),POINTER  ::  fc  
      REAL(ppm_kind_double), DIMENSION(:,:,:),POINTER  ::  err  
      !lets save the boundary condition.index:face of the subdomain(1:6)
      TYPE(bc_value_3d_sca_d), DIMENSION(:), POINTER   ::  bcvalue   
   END TYPE mg_field_3d_sca_d
#undef  __KIND


#define __KIND == __SINGLE_PRECISION
   !1st index: subdomain,2nd index : multigrid level
   TYPE(mg_field_3d_sca_s), DIMENSION(:,:),              POINTER  ::   mgfield_3d_sca_s 
#undef __KIND
#define __KIND == __DOUBLE_PRECISION
   !1st index: subdomain,2nd index : multigrid level
   TYPE(mg_field_3d_sca_d), DIMENSION(:,:),              POINTER  ::   mgfield_3d_sca_d 
#undef __KIND
#undef __MESH_DIM

#undef __DIM

#define __DIM __VFIELD
#define __MESH_DIM __2D
#define __KIND __SINGLE_PRECISION
   TYPE bc_value_2d_vec_s
      ! 1st index mesh position locally
      REAL(ppm_kind_single), DIMENSION(:),POINTER   ::  pbcvalue
   END TYPE bc_value_2d_vec_s
#undef  __KIND
#define __KIND == __DOUBLE_PRECISION
   TYPE bc_value_2d_vec_d
      !1st index mesh position locally
      REAL(ppm_kind_single), DIMENSION(:),POINTER   ::   pbcvalue
   END TYPE bc_value_2d_vec_d
#undef __KIND

#define __KIND == __SINGLE_PRECISION
   !----------------------------------------------------------------------------
   ! Our multigrid field with all its necessary components (Take a look at the 
   ! theory)
   !----------------------------------------------------------------------------
   TYPE mg_field_2d_vec_s
      !function corrections, error restrictions, errors
      !1st index component 2nd and 3rd index: mesh position(local)
      REAL(ppm_kind_single), DIMENSION(:,:,:),POINTER  ::  uc 
      REAL(ppm_kind_single), DIMENSION(:,:,:),POINTER  ::  fc 
      REAL(ppm_kind_single), DIMENSION(:,:,:),POINTER  ::  err  

      !lets save the boundary condition.index:component,face of the subdomain(1:4)
      TYPE(bc_value_2d_vec_s), DIMENSION(:,:), POINTER   ::  bcvalue   
   END TYPE mg_field_2d_vec_s
#undef  __KIND


#define __KIND  == __DOUBLE_PRECISION  
   TYPE mg_field_2d_vec_d
      !function corrections, error restrictions, errors,
      !1st index: component 2nd:3rd index: mesh position(local)
      REAL(ppm_kind_double), DIMENSION(:,:,:),POINTER  ::  uc  
      REAL(ppm_kind_double), DIMENSION(:,:,:),POINTER  ::  fc  
      REAL(ppm_kind_double), DIMENSION(:,:,:),POINTER  ::  err  
      !lets save the boundary condition.index:component,face of the subdomain(1:4)
      TYPE(bc_value_2d_vec_d), DIMENSION(:), POINTER   ::  bcvalue   
   END TYPE mg_field_2d_vec_d
#undef  __KIND


#define __KIND == __SINGLE_PRECISION
   !1st index: subdomain,2nd index : multigrid level
   TYPE(mg_field_2d_vec_s), DIMENSION(:,:),              POINTER  ::   mgfield_2d_vec_s 
#undef __KIND
#define __KIND == __DOUBLE_PRECISION
   !1st index: subdomain,2nd index : multigrid level
   TYPE(mg_field_2d_vec_d), DIMENSION(:,:),              POINTER  ::   mgfield_2d_vec_d 
#undef __KIND
#undef __MESH_DIM

#define __MESH_DIM __3D
#define __KIND __SINGLE_PRECISION
   TYPE bc_value_3d_vec_s
      ! 1st-2nd index mesh position locally
      REAL(ppm_kind_single), DIMENSION(:,:,:),POINTER   ::  pbcvalue
   END TYPE bc_value_3d_vec_s
#undef  __KIND
#define __KIND == __DOUBLE_PRECISION
   TYPE bc_value_3d_vec_d
      !1st-2nd index mesh position locally
      REAL(ppm_kind_single), DIMENSION(:,:,:),POINTER   ::   pbcvalue
   END TYPE bc_value_3d_vec_d
#undef __KIND



#define __KIND == __SINGLE_PRECISION
   !----------------------------------------------------------------------------
   ! Our multigrid field with all its necessary components (Take a look at the 
   ! theory)
   !----------------------------------------------------------------------------
   TYPE mg_field_3d_vec_s
      !function corrections, error restrictions, errors
      !1st index: component 2nd-4th index: mesh position(local)
      REAL(ppm_kind_single), DIMENSION(:,:,:,:),POINTER  ::  uc 
      REAL(ppm_kind_single), DIMENSION(:,:,:,:),POINTER  ::  fc 
      REAL(ppm_kind_single), DIMENSION(:,:,:,:),POINTER  ::  err  

      !lets save the boundary condition.index:component,face of the subdomain(1:6)
      TYPE(bc_value_3d_vec_s), DIMENSION(:), POINTER   ::  bcvalue   
   END TYPE mg_field_3d_vec_s
#undef  __KIND

#define __KIND  == __DOUBLE_PRECISION  

   TYPE mg_field_3d_vec_d
      !function corrections, error restrictions, errors,
      !1st index component,2nd:4th index: mesh position(local)
      REAL(ppm_kind_double), DIMENSION(:,:,:,:),POINTER  ::  uc  
      REAL(ppm_kind_double), DIMENSION(:,:,:,:),POINTER  ::  fc  
      REAL(ppm_kind_double), DIMENSION(:,:,:,:),POINTER  ::  err  
      !lets save the boundary condition.index:face of the subdomain(1:6)
      TYPE(bc_value_3d_vec_d), DIMENSION(:), POINTER   ::  bcvalue   
   END TYPE mg_field_3d_vec_d
#undef  __KIND


#define __KIND == __SINGLE_PRECISION
   !1st index: subdomain,2nd index : multigrid level
   TYPE(mg_field_3d_vec_s), DIMENSION(:,:),              POINTER  ::   mgfield_3d_vec_s 
#undef __KIND
#define __KIND == __DOUBLE_PRECISION
   !1st index: subdomain,2nd index : multigrid level
   TYPE(mg_field_3d_vec_d), DIMENSION(:,:),              POINTER  ::   mgfield_3d_vec_d 
#undef __KIND
#undef __MESH_DIM

#undef __DIM
   !----------------------------------------------------------------------------
   !Starting  index for the iteration through the mesh points. 
   !----------------------------------------------------------------------------
   INTEGER,  DIMENSION(:,:,:), POINTER  :: start
   !----------------------------------------------------------------------------
   !Stopping index for the iteration through the mesh points.
   !----------------------------------------------------------------------------
   INTEGER,  DIMENSION(:,:,:), POINTER  :: istop 
   !----------------------------------------------------------------------------
   !Factor for coarsening the mesh
   !----------------------------------------------------------------------------
   INTEGER,  DIMENSION(:),POINTER          :: factor
   !----------------------------------------------------------------------------
   !Array with internal meshids
   !----------------------------------------------------------------------------
   INTEGER,  DIMENSION(:),POINTER          :: meshid_g
   !----------------------------------------------------------------------------
   !Array with external mesh_ids
   !----------------------------------------------------------------------------
   INTEGER,  DIMENSION(:),POINTER          :: mesh_id_g
   !----------------------------------------------------------------------------
   !Size of the ghostlayer. It is 1 for the multigrid since we do
   !for the time being second order finite differences
   !----------------------------------------------------------------------------
   INTEGER,  DIMENSION(:),POINTER       :: ghostsize
   !----------------------------------------------------------------------------
   !BOUNDARY CONDITIONS of the computational domain.1st index:sub,2nd,face
   !----------------------------------------------------------------------------
#define __DIM == __SFIELD
   INTEGER,  DIMENSION(:,:),POINTER     :: bcdef_sca
#undef __DIM
#define __DIM == __VFIELD
   INTEGER,  DIMENSION(:,:,:),POINTER     :: bcdef_vec
#undef __DIM

   !----------------------------------------------------------------------------
   !Is the face of the cell at the boundary? Yes or no?1st index face,2nd:isub
   !----------------------------------------------------------------------------
   LOGICAL,  DIMENSION(:,:), POINTER  :: lboundary
   !----------------------------------------------------------------------------
   !V_CYCLE OR W_CYCLE AND TO PRINT OR NOT TO PRINT
   !----------------------------------------------------------------------------
   LOGICAL                             :: w_cycle
   LOGICAL                             :: l_print
   !----------------------------------------------------------------------------
   !ARE ALL THE BOUNDARIES PERIODIC
   !----------------------------------------------------------------------------
   LOGICAL                             :: lperiodic
   !----------------------------------------------------------------------------
   !Order of the mg
   !----------------------------------------------------------------------------
   INTEGER                             :: order
   !----------------------------------------------------------------------------
   !number of levels (theoretical value)
   !----------------------------------------------------------------------------
   INTEGER                              :: maxlev
   !----------------------------------------------------------------------------
   !number of subs
   !----------------------------------------------------------------------------
   INTEGER                              :: nsubs
   !----------------------------------------------------------------------------
   !smoother              
   !----------------------------------------------------------------------------
   INTEGER                              :: ismoother
   !----------------------------------------------------------------------------
   !number of dimensions in the problem(if scalar fields=> vecdim=1)
   !----------------------------------------------------------------------------
   INTEGER                              :: vecdim

   !----------------------------------------------------------------------------
   !Array with the maximum number of mesh points on each processor
   !Due to the load ballancing the waste of memory (if existed) is 
   !minimal !!
   !----------------------------------------------------------------------------
   INTEGER,DIMENSION(:,:),POINTER :: max_node


#define __KIND __SINGLE_PRECISION 
   REAL(ppm_kind_single)                        :: rdx2_s,rdy2_s,rdz2_s
   REAL(ppm_kind_single)                        :: dx_s,dy_s,dz_s
   REAL(ppm_kind_single)                        :: EPSU_s
   REAL(ppm_kind_single)                        :: omega_s
#undef __KIND

#define __KIND __DOUBLE_PRECISION 
   REAL(ppm_kind_double)                        :: rdx2_d,rdy2_d,rdz2_d
   REAL(ppm_kind_double)                        :: dx_d,dy_d,dz_d
   REAL(ppm_kind_double)                        :: EPSU_d
   REAL(ppm_kind_double)                        :: omega_d
#undef __KIND

 END MODULE ppm_module_data_mg


