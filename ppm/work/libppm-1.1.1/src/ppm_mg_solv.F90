#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

       !------------------------------------------------------------------------
       !  Subroutine   :                  ppm_mg_solv 
       !------------------------------------------------------------------------
       !
       !  Input        :    itera      (I)  :  initial smoothing sweeps 
       !                                       in the finest level.
       !              
       !                    iterf      (I)  :  final smoothing sweeps
       !                                       in the finest level
       !
       !                    iter1      (I)  :  AFTER EACH RESTRICTION   
       !                                       SMOOTHING SWEEPS TAKE PLACE
       !                                       IMPORTANT PARAMETER
       !
       !                    iter2      (I)  :  AFTER EACH PROLONGATION
       !                                       SMOOTHING SWEEPS TAKE PLACE
       !                                        
       !
       !  Input/Output :     u         (F)  :  THE FIELD OF THE SOLUTION
       !                                       WITH GHOST VALUES!!
       !                     f         (F)  :  THE FIELD OF THE RHS (NO GHOST
       !                                         VALUES)
       !  Output       :    info       (I)
       !
       !  Purpose      : 
       !
       !
       !  References   :
       !
       !  Revisions    :
       !------------------------------------------------------------------------
       !  $Log: ppm_mg_solv.f,v $
       !  Revision 1.16  2006/07/21 11:30:54  kotsalie
       !  FRIDAY
       !
       !  Revision 1.14  2005/12/08 12:44:46  kotsalie
       !  commiting dirichlet
       !
       !  Revision 1.13  2005/05/30 13:03:22  kotsalie
       !  UPDATED FOR SERIAL VERSION WITHOUT MPI
       !
       !  Revision 1.12  2005/03/14 13:24:03  kotsalie
       !  COMMITED THE VECTOR CASE. IT IS FOR LDA=3
       !
       !  Revision 1.11  2005/01/04 09:48:21  kotsalie
       !  ghostsize=2 scalar case
       !
       !  Revision 1.10  2004/11/05 15:18:35  kotsalie
       !  Made independent the initial and final smoothing steps
       !
       !  Revision 1.9  2004/10/13 16:02:03  kotsalie
       !  Maximum residual between processors is communicated
       !
       !  Revision 1.8  2004/09/30 14:26:24  kotsalie
       !  *** empty log message ***
       !
       !  Revision 1.7  2004/09/29 10:47:36  kotsalie
       !  The user can now print the residual. THis should serve for him
       !  as a istopping criterium
       !
       !  Revision 1.6  2004/09/23 13:50:54  kotsalie
       !  Changed IF (w_cycle) to IF(.FALSE.) Now the recusrion goes up to level 2.
       !
       !  Revision 1.5  2004/09/23 12:41:16  kotsalie
       !  MG new version
       !
       !------------------------------------------------------------------------
       !  Parallel Particle Mesh Library (PPM)
       !  Institute of Computational Science
       !  ETH Zentrum, Hirschengraben 84
       !  CH-8092 Zurich, Switzerland
       !------------------------------------------------------------------------

#if   __DIM   == __SFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_solv_2d_sca_s(u,f,itera,iterf,iter1,iter2,Eu,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_solv_2d_sca_d(u,f,itera,iterf,iter1,iter2,Eu,info)
#endif 
#elif __MESH_DIM   == __3D
#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_solv_3d_sca_s(u,f,itera,iterf,iter1,iter2,Eu,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_solv_3d_sca_d(u,f,itera,iterf,iter1,iter2,Eu,info)
#endif 
#endif 
#elif __DIM == __VFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_solv_2d_vec_s(u,f,lda,itera,iterf,iter1,iter2,Eu,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_solv_2d_vec_d(u,f,lda,itera,iterf,iter1,iter2,Eu,info)
#endif
#elif __MESH_DIM   == __3D
#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_solv_3d_vec_s(u,f,lda,itera,iterf,iter1,iter2,Eu,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_solv_3d_vec_d(u,f,lda,itera,iterf,iter1,iter2,Eu,info)
#endif
#endif
#endif

#include "ppm_define.h"

          !---------------------------------------------------------------------
         !  Modules 
         !----------------------------------------------------------------------
         USE ppm_module_data
         USE ppm_module_data_mg
         USE ppm_module_data_mesh
         USE ppm_module_substart
         USE ppm_module_substop
         USE ppm_module_error
         USE ppm_module_alloc
         USE ppm_module_map_field_ghost
         USE ppm_module_mg_core
         USE ppm_module_mg_res
         USE ppm_module_mg_prolong         
         USE ppm_module_mg_smooth         
         USE ppm_module_write

         IMPLICIT NONE

#ifdef HAVE_MPI
       INCLUDE  'mpif.h'
#else
#include "fakempi.h"
#endif

#if    __KIND == __SINGLE_PRECISION
         INTEGER, PARAMETER :: MK = ppm_kind_single
#else
         INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
         !----------------------------------------------------------------------
         !  Arguments (for u and f index: local mesh locations and isub) 
         !----------------------------------------------------------------------
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:),POINTER     ::  u
         REAL(MK),DIMENSION(:,:,:),POINTER     ::  f
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:),POINTER   ::  u
         REAL(MK),DIMENSION(:,:,:,:),POINTER   ::  f
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:,:),POINTER     ::  u
         REAL(MK),DIMENSION(:,:,:,:),POINTER     ::  f
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:,:),POINTER   ::  u
         REAL(MK),DIMENSION(:,:,:,:,:),POINTER   ::  f
#endif
#endif
#if __DIM == __VFIELD
         INTEGER,INTENT(IN)                      :: lda           
#endif
         INTEGER,                   INTENT(IN)   ::  itera,iterf,iter1,iter2
         REAL(MK),                  INTENT(OUT)  ::  Eu  
         INTEGER,                   INTENT(INOUT)   ::  info
         !----------------------------------------------------------------------
         !  Local variables 
         !----------------------------------------------------------------------
         REAL(MK)                             :: t0
         REAL(MK)                             :: E,res
         INTEGER                              :: iface,count,k
         INTEGER                              :: ix,iy  
         CHARACTER(LEN=256)                   :: cbuf
         INTEGER                              :: mlev,color,it
         INTEGER                              :: ncalls=0
         REAL(MK)                             :: c1,c2,c3,c4  
         INTEGER                              :: isub,i,j
         REAL(MK)                             :: x,y
         REAL(MK)                             :: gEu 
         INTEGER                              :: MPI_PREC


#if __MESH_DIM == __3D
         REAL(MK)                             :: c5,dz,rdz2
         INTEGER,DIMENSION(4)                 :: ldl4,ldu4
         INTEGER,DIMENSION(5)                 :: ldl5,ldu5
#endif
         INTEGER                              :: ilda
         REAL(MK)                             :: rdx2,rdy2
         REAL(MK)                             :: dx,dy
#if __MESH_DIM == __2D
         INTEGER,DIMENSION(3)                 :: ldl3,ldu3
         INTEGER,DIMENSION(4)                 :: ldl4,ldu4
#endif
         INTEGER                              :: topoid,iopt,idom

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         TYPE(mg_field_2d_sca_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_2d_sca_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
         TYPE(mg_field_3d_sca_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_3d_sca_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         TYPE(mg_field_2d_vec_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_2d_vec_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
         TYPE(mg_field_3d_vec_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
         TYPE(mg_field_3d_vec_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#endif
#endif


#if __DIM == __SFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:),POINTER :: uc_dummy
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:),POINTER :: uc_dummy
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:,:),POINTER :: uc_dummy
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:,:),POINTER :: uc_dummy
#endif
#endif


#if __DIM == __SFIELD
#if __MESH_DIM == __2D
     REAL(MK),DIMENSION(:,:),POINTER :: tuc
#elif __MESH_DIM == __3D
     REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
     REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
#elif __MESH_DIM == __3D
     REAL(MK),DIMENSION(:,:,:,:),POINTER :: tuc
#endif
#endif
         !----------------------------------------------------------------------
         !  Externals 
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         !  Initialize 
         !----------------------------------------------------------------------

         CALL substart('ppm_mg_solv',t0,info)


#ifdef USE_MPI
         IF (ppm_kind.EQ.ppm_kind_single) THEN
            MPI_PREC = MPI_REAL
         ELSE
            MPI_PREC = MPI_DOUBLE_PRECISION
         ENDIF
#endif

         !----------------------------------------------------------------------
         !  Check arguments
         !----------------------------------------------------------------------
         IF (ppm_debug .GT. 0) THEN
#if __DIM == __SFIELD        
#if __MESH_DIM == __2D        
            IF (SIZE(u,3) .LT. nsubs) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution exist on nsubs subdomains',__LINE__,info)        
               GOTO 9999
            ENDIF
            topoid=ppm_field_topoid
            DO i=1,nsubs
               idom=ppm_isublist(i,topoid)
               IF (SIZE(u(:,:,i),1).LT.ppm_cart_mesh(meshid_g(1),  &
      &           topoid)%nnodes(1,idom)+2*ghostsize(1)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution mess with mesh points in x-dir!',__LINE__,info)  
                  GOTO 9999    
               ENDIF
               IF (SIZE(u(:,:,i),2).LT.ppm_cart_mesh(meshid_g(1),  &
      &           topoid)%nnodes(2,idom)+2*ghostsize(2)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution mess with mesh points in y-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
            ENDDO
            IF (SIZE(f,3) .LT. nsubs) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'rhs exist on nsubs subdomains!',__LINE__,info)  
               GOTO 9999
            ENDIF
            topoid=ppm_field_topoid
            DO i=1,nsubs
               idom=ppm_isublist(i,topoid)
               IF (SIZE(f(:,:,i),1).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                  topoid)%nnodes(1,idom)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'rhs mess with mesh points in x-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
               IF (SIZE(f(:,:,i),2).LT.ppm_cart_mesh(meshid_g(1),  &
      &           topoid)%nnodes(2,idom)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'rhs mess with mesh points in y-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
            ENDDO
#elif __MESH_DIM == __3D
            IF (SIZE(u,4) .LT. nsubs) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution exist on nsubs subdomains!',__LINE__,info)        
               GOTO 9999
            ENDIF
            topoid=ppm_field_topoid
            DO i=1,nsubs
               idom=ppm_isublist(i,topoid)
               IF (SIZE(u(:,:,:,i),1).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                   topoid)%nnodes(1,idom)+2*ghostsize(1)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution mess with mesh points in x-dir!',__LINE__,info)  
                  GOTO 9999    
               ENDIF
               IF (SIZE(u(:,:,:,i),2).LT.ppm_cart_mesh(meshid_g(1),  &
      &           topoid)%nnodes(2,idom)+2*ghostsize(2)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution mess with mesh points in y-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
               IF (SIZE(u(:,:,:,i),3).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                   topoid)%nnodes(3,idom)+2*ghostsize(3)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution mess with mesh points in z-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
           ENDDO
            IF (SIZE(f,4) .LT. nsubs) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'rhs exist on nsubs subdomains!',__LINE__,info)  
               GOTO 9999
            ENDIF
            topoid=ppm_field_topoid
            DO i=1,nsubs
               idom=ppm_isublist(i,topoid)
               IF (SIZE(f(:,:,:,i),1).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                    topoid)%nnodes(1,idom)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &                   'rhs mess with mesh points in x-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
               IF (SIZE(f(:,:,:,i),2).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                          topoid)%nnodes(2,idom)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &            'rhs mess with mesh points in y-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
               IF (SIZE(f(:,:,:,i),3).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                         topoid)%nnodes(3,idom)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'rhs mess with mesh points in z-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
            ENDDO
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D        
            IF (SIZE(u,4) .LT. nsubs) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution exist on nsubs subdomains',__LINE__,info)        
               GOTO 9999
            ENDIF
            topoid=ppm_field_topoid
            DO i=1,nsubs
               idom=ppm_isublist(i,topoid)
               IF (SIZE(u(:,:,:,i),2).LT.ppm_cart_mesh(meshid_g(1),  &
      &           topoid)%nnodes(1,idom)+2) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution mess with mesh points in x-dir!',__LINE__,info)  
                  GOTO 9999    
               ENDIF
               IF (SIZE(u(:,:,:,i),3).LT.ppm_cart_mesh(meshid_g(1),  &
      &           topoid)%nnodes(2,idom)+2) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution mess with mesh points in y-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
            ENDDO
            IF (SIZE(f,4) .LT. nsubs) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'rhs exist on nsubs subdomains!',__LINE__,info)  
               GOTO 9999
            ENDIF
            topoid=ppm_field_topoid
            DO i=1,nsubs
               idom=ppm_isublist(i,topoid)
               IF (SIZE(f(:,:,:,i),2).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                  topoid)%nnodes(1,idom)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'rhs mess with mesh points in x-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
               IF (SIZE(f(:,:,:,i),3).LT.ppm_cart_mesh(meshid_g(1),  &
      &           topoid)%nnodes(2,idom)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'rhs mess with mesh points in y-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
            ENDDO
#elif __MESH_DIM == __3D
            IF (SIZE(u,5) .LT. nsubs) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution exist on nsubs subdomains!',__LINE__,info)        
               GOTO 9999
            ENDIF
            topoid=ppm_field_topoid
            DO i=1,nsubs
               idom=ppm_isublist(i,topoid)
               IF (SIZE(u(:,:,:,:,i),2).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                   topoid)%nnodes(1,idom)+2*ghostsize(1)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution mess with mesh points in x-dir!',__LINE__,info)  
                  GOTO 9999    
               ENDIF
               IF (SIZE(u(:,:,:,:,i),3).LT.ppm_cart_mesh(meshid_g(1),  &
      &           topoid)%nnodes(2,idom)+2*ghostsize(2)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution mess with mesh points in y-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
               IF (SIZE(u(:,:,:,:,i),4).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                   topoid)%nnodes(3,idom)+2*ghostsize(3)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'solution mess with mesh points in z-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
           ENDDO
            IF (SIZE(f,5) .LT. nsubs) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'rhs exist on nsubs subdomains!',__LINE__,info)  
               GOTO 9999
            ENDIF
            topoid=ppm_field_topoid
            DO i=1,nsubs
               idom=ppm_isublist(i,topoid)
               IF (SIZE(f(:,:,:,:,i),2).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                    topoid)%nnodes(1,idom)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &                   'rhs mess with mesh points in x-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
               IF (SIZE(f(:,:,:,:,i),3).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                          topoid)%nnodes(2,idom)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &            'rhs mess with mesh points in y-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
               IF (SIZE(f(:,:,:,:,i),4).LT.ppm_cart_mesh(meshid_g(1),  &
      &                                         topoid)%nnodes(3,idom)) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_solv',  &
      &             'rhs mess with mesh points in z-dir!',__LINE__,info)
                  GOTO 9999
               ENDIF
            ENDDO
#endif
#endif
         ENDIF

         !----------------------------------------------------------------------
         !Definition of necessary variables and allocation of arrays
         !----------------------------------------------------------------------
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
#if __DIM == __SFIELD
         mgfield=>mgfield_2d_sca_s
#elif __DIM == __VFIELD
         mgfield=>mgfield_2d_vec_s
#endif
         rdx2=rdx2_s
         rdy2=rdy2_s
         dx=dx_s
         dy=dy_s
#elif __KIND == __DOUBLE_PRECISION
#if __DIM == __SFIELD
         mgfield=>mgfield_2d_sca_d
#elif __DIM == __VFIELD
         mgfield=>mgfield_2d_vec_d
#endif
         rdx2=rdx2_d
         rdy2=rdy2_d
         dx=dx_d
         dy=dy_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
#if __DIM == __SFIELD
         mgfield=>mgfield_3d_sca_s
#elif __DIM == __VFIELD 
         mgfield=>mgfield_3d_vec_s
#endif
         rdx2=rdx2_s
         rdy2=rdy2_s
         rdz2=rdz2_s
         dx=dx_s
         dy=dy_s
         dz=dz_s
#elif __KIND == __DOUBLE_PRECISION
#if __DIM == __SFIELD
         mgfield=>mgfield_3d_sca_d
#elif __DIM == __VFIELD
         mgfield=>mgfield_3d_vec_d
#endif
         rdx2=rdx2_d
         rdy2=rdy2_d
         rdz2=rdz2_d
         dx=dx_d
         dy=dy_d
         dz=dz_d
#endif
#endif

         topoid=ppm_field_topoid

      ncalls=ncalls+1
      IF (ncalls.EQ.1) THEN

         DO i=1,maxlev

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
              iopt = ppm_param_alloc_fit
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1
              ldu3(1) = max_node(1,i)+ghostsize(1)
              ldu3(2) = max_node(2,i)+ghostsize(2)
              ldu3(3) = nsubs
              CALL ppm_alloc(uc_dummy,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'MGsolv',    &
      &                       'uc_dummy',__LINE__,info)
                GOTO 9999
               ENDIF
              uc_dummy(:,:,:)=0.0_MK

#elif __MESH_DIM ==__3D
              iopt = ppm_param_alloc_fit
              ldl4(1) = 1-ghostsize(1)
              ldl4(2) = 1-ghostsize(2)
              ldl4(3) = 1-ghostsize(3)
              ldl4(4) = 1
              ldu4(1) = max_node(1,i)+ghostsize(1)
              ldu4(2) = max_node(2,i)+ghostsize(2)
              ldu4(3) = max_node(3,i)+ghostsize(3)
              ldu4(4) = nsubs
              CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'MGsolv',    &
       &                       'uc_dummy',__LINE__,info)
                GOTO 9999
               ENDIF

              uc_dummy(:,:,:,:)=0.0_MK
#endif

            CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(i),&
      &                             ghostsize,ppm_param_map_init,info) 




#if __MESH_DIM == __2D
              iopt = ppm_param_dealloc
              ldl3(1) = 1-ghostsize(1)
              ldl3(2) = 1-ghostsize(2)
              ldl3(3) = 1
              ldu3(1) = max_node(1,i)+ghostsize(1)
              ldu3(2) = max_node(2,i)+ghostsize(2)
              ldu3(3) = nsubs
              CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'MGsolv',    &
      &                       'uc_dummy',__LINE__,info)
                GOTO 9999
               ENDIF

#elif __MESH_DIM ==__3D
              iopt = ppm_param_dealloc
              ldl4(1) = 1-ghostsize(1)
              ldl4(2) = 1-ghostsize(2)
              ldl4(3) = 1-ghostsize(3)
              ldl4(4) = 1
              ldu4(1) = max_node(1,i)+ghostsize(1)
              ldu4(2) = max_node(2,i)+ghostsize(2)
              ldu4(3) = max_node(3,i)+ghostsize(3)
              ldu4(4) = nsubs
              CALL ppm_alloc(uc_dummy,ldl5,ldu5,iopt,info)
              IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'MGsolv',    &
       &                       'uc_dummy',__LINE__,info)
                GOTO 9999
               ENDIF

#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
              iopt = ppm_param_alloc_fit
              ldl4(1) = 1
              ldl4(2) = 1-ghostsize(1)
              ldl4(3) = 1-ghostsize(2)
              ldl4(4) = 1
              ldu4(1) = vecdim
              ldu4(2) = max_node(1,i)+ghostsize(1)
              ldu4(3) = max_node(2,i)+ghostsize(2)
              ldu4(4) = nsubs
              CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'MGsolv',    &
      &                       'uc_dummy',__LINE__,info)
                GOTO 9999
               ENDIF
              uc_dummy(:,:,:,:)=0.0_MK

#elif __MESH_DIM ==__3D
              iopt = ppm_param_alloc_fit
              ldl5(1) = 1
              ldl5(2) = 1-ghostsize(1)
              ldl5(3) = 1-ghostsize(2)
              ldl5(4) = 1-ghostsize(3)
              ldl5(5) = 1
              ldu5(1) = vecdim
              ldu5(2) = max_node(1,i)+ghostsize(1)
              ldu5(3) = max_node(2,i)+ghostsize(2)
              ldu5(4) = max_node(3,i)+ghostsize(3)
              ldu5(5) = nsubs
              CALL ppm_alloc(uc_dummy,ldl5,ldu5,iopt,info)
              IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'MGsolv',    &
       &                       'uc_dummy',__LINE__,info)
                GOTO 9999
               ENDIF

              uc_dummy(:,:,:,:,:)=0.0_MK
#endif

            CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(i),&
      &                             ghostsize,ppm_param_map_init,info) 




#if __MESH_DIM == __2D
              iopt = ppm_param_dealloc
              ldl4(1) = 1-ghostsize(1)
              ldl4(1) = 1-ghostsize(2)
              ldl4(1) = 1
              ldu4(1) = max_node(1,i)+ghostsize(1)
              ldu4(2) = max_node(2,i)+ghostsize(2)
              ldu4(3) = nsubs
              CALL ppm_alloc(uc_dummy,ldl3,ldu3,iopt,info)
              IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'MGsolv',    &
      &                       'uc_dummy',__LINE__,info)
                GOTO 9999
               ENDIF

#elif __MESH_DIM ==__3D
              iopt = ppm_param_dealloc
              ldl5(1) = 1-ghostsize(1)
              ldl5(2) = 1-ghostsize(2)
              ldl5(3) = 1-ghostsize(3)
              ldl5(4) = 1
              ldu5(1) = max_node(1,i)+ghostsize(1)
              ldu5(2) = max_node(2,i)+ghostsize(2)
              ldu5(3) = max_node(3,i)+ghostsize(3)
              ldu5(4) = nsubs
              CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
              IF (info .NE. 0) THEN
                info = ppm_error_fatal
                CALL ppm_error(ppm_err_alloc,'MGsolv',    &
       &                       'uc_dummy',__LINE__,info)
                GOTO 9999
               ENDIF

#endif
#endif



         ENDDO

         ncalls=ncalls+1

       ENDIF

         !----------------------------------------------------------------------
         !DO n1 initial sweeps in the finest mesh  with a GS-solver to get the 
         !initial solution 
         !----------------------------------------------------------------------
#if __DIM == __SFIELD
#if __MESH_DIM == __2D


         c1 = 1.0_MK/(2.0_MK*(rdx2+rdy2))  
         c2 = rdx2
         c3 = rdy2     
         c4 = 2.0_MK*c2+2.0_MK*c3
         count = 0


         CALL ppm_mg_smooth_sca(u,f,itera,1,c1,c2,c3,info) 

         !----------------------------------------------------------------------
         ! Compute residual
         !----------------------------------------------------------------------

         CALL ppm_mg_res_sca(u,f,c1,c2,c3,c4,E,info)

#ifdef USE_MPI
         CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
         E=gEu
#endif

         IF (info .NE. 0) THEN 
          GOTO 9999
         ENDIF 

       IF (l_print) THEN 
         WRITE(cbuf,*) 'Eu:',E
         CALL PPM_WRITE(ppm_rank,'mg_solv',cbuf,info)
       ENDIF

          !---------------------------------------------------------------------
         !Initiation of the function correction. (We start on purpose with lev=2)
         !----------------------------------------------------------------------
         DO mlev=2,maxlev
            DO isub=1,nsubs
               tuc=>mgfield(isub,mlev)%uc
               DO j=start(2,isub,mlev),istop(2,isub,mlev)
                  DO i=start(1,isub,mlev),istop(1,isub,mlev)
                        tuc(i,j)=0.0_MK
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !----------------------------------------------------------------------
         !CALL THE MULTIGRID TO DO NICE STUFF TO OUR FUNCTION
         !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
         CALL ppm_mg_core_2d_sca_s(2,iter1,iter2,info)  
#elif __KIND == __DOUBLE_PRECISION
         CALL ppm_mg_core_2d_sca_d(2,iter1,iter2,info)  
#endif   

         !----------------------------------------------------------------------
         !PROLONG the solution to the finest grid
         !----------------------------------------------------------------------

#if __KIND == __SINGLE_PRECISION
          CALL ppm_mg_prolong_2d_sca_s(1,info)
#elif __KIND == __DOUBLE_PRECISION
          CALL ppm_mg_prolong_2d_sca_d(1,info)
#endif

         !----------------------------------------------------------------------
         !UPDATE THE FUNCTION
         !----------------------------------------------------------------------
         DO isub=1,nsubs
            tuc=>mgfield(isub,mlev)%uc
            DO j=start(2,isub,1),istop(2,isub,1)   
               DO i=start(1,isub,1),istop(1,isub,1)
                     u(i,j,isub)=tuc(i,j) 
               ENDDO
            ENDDO
         ENDDO

         !----------------------------------------------------------------------
         !DO the final sweeps
           !--------------------------------------------------------------------


         CALL ppm_mg_smooth_sca(u,f,iterf,1,c1,c2,c3,info) 
         CALL ppm_mg_res_sca(u,f,c1,c2,c3,c4,E,info)
#ifdef USE_MPI
         CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
         Eu=gEu
#else
         Eu=E
#endif


#elif __MESH_DIM == __3D



         c1 = 1.0_MK/(2.0_MK*(rdx2+rdy2+rdz2))
         c2 = rdx2
         c3 = rdy2
         c4 = rdz2 
         c5 = 2.0_MK*c2+2.0_MK*c3+2.0_MK*c4

         CALL ppm_mg_smooth_sca(u,f,itera,1,c1,c2,c3,c4,info) 

         !----------------------------------------------------------------------
         ! Compute residual
         !----------------------------------------------------------------------


         CALL ppm_mg_res_sca(u,f,c1,c2,c3,c4,c5,E,info)

#ifdef USE_MPI
         CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
         E=gEu
#endif

       IF (l_print) THEN 
         WRITE(cbuf,*) 'Eu:',E
         CALL PPM_WRITE(ppm_rank,'mg_solv',cbuf,info)
       ENDIF
          !---------------------------------------------------------------------
         !Initiation of the function correction. (We start on purpose with lev=2)
         !----------------------------------------------------------------------
         DO mlev=2,maxlev
            DO isub=1,nsubs
               tuc=>mgfield(isub,mlev)%uc
               DO k=start(3,isub,mlev),istop(3,isub,mlev) 
                  DO j=start(2,isub,mlev),istop(2,isub,mlev)
                     DO i=start(1,isub,mlev),istop(1,isub,mlev)
                           tuc(i,j,k)=0.0_MK
                     ENDDO
                 ENDDO
               ENDDO
            ENDDO
         ENDDO
         !----------------------------------------------------------------------
         !CALL THE MULTIGRID TO DO NICE STUFF TO OUR FUNCTION
         !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
         CALL ppm_mg_core_3d_sca_s(2,iter1,iter2,info)
#elif __KIND == __DOUBLE_PRECISION
         CALL ppm_mg_core_3d_sca_d(2,iter1,iter2,info)
#endif
         !----------------------------------------------------------------------
         !PROLONG the solution to the finest grid
         !----------------------------------------------------------------------

#if __KIND == __SINGLE_PRECISION
         CALL ppm_mg_prolong_3d_sca_s(1,info)
#elif __KIND == __DOUBLE_PRECISION
         CALL ppm_mg_prolong_3d_sca_d(1,info)
#endif

         !----------------------------------------------------------------------
         !UPDATE THE FUNCTION
         !----------------------------------------------------------------------
         DO isub=1,nsubs
               tuc=>mgfield(isub,mlev)%uc
            DO k=start(3,isub,1),istop(3,isub,1)
               DO j=start(2,isub,1),istop(2,isub,1)
                  DO i=start(1,isub,1),istop(1,isub,1)
                        u(i,j,k,isub)=tuc(i,j,k)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         !----------------------------------------------------------------------
         !DO the final sweeps
           !--------------------------------------------------------------------
         CALL ppm_mg_smooth_sca(u,f,iterf,1,c1,c2,c3,c4,info) 
         CALL ppm_mg_res_sca(u,f,c1,c2,c3,c4,c5,E,info)

#ifdef USE_MPI
         CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
         Eu=gEu
#else
         Eu=E
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D


         c1 = 1.0_MK/(2.0_MK*(rdx2+rdy2))  
         c2 = rdx2
         c3 = rdy2     
         c4 = 2.0_MK*c2+2.0_MK*c3
         count = 0


         CALL ppm_mg_smooth_vec(u,f,itera,1,c1,c2,c3,info) 


         !----------------------------------------------------------------------
         ! Compute residual
         !----------------------------------------------------------------------

         CALL ppm_mg_res_vec(u,f,c1,c2,c3,c4,E,info)   

#ifdef USE_MPI
         CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
         E=gEu
#endif


         IF (l_print) THEN 
          WRITE(cbuf,*) 'Eu:',E
          CALL PPM_WRITE(ppm_rank,'mg_solv',cbuf,info)
         ENDIF

          !---------------------------------------------------------------------
         !Initiation of the function correction. (We start on purpose with lev=2)
         !----------------------------------------------------------------------
         DO mlev=2,maxlev
            DO isub=1,nsubs
               tuc=>mgfield(isub,mlev)%uc
               DO j=start(2,isub,mlev),istop(2,isub,mlev)
                  DO i=start(1,isub,mlev),istop(1,isub,mlev)
                   DO ilda=1,vecdim
                        tuc(ilda,i,j)=0.0_MK
                   ENDDO 
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !----------------------------------------------------------------------
         !CALL THE MULTIGRID TO DO NICE STUFF TO OUR FUNCTION
         !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
         CALL ppm_mg_core_2d_vec_s(2,iter1,iter2,info)  
#elif __KIND == __DOUBLE_PRECISION
         CALL ppm_mg_core_2d_vec_d(2,iter1,iter2,info)  
#endif   

         !----------------------------------------------------------------------
         !PROLONG the solution to the finest grid
         !----------------------------------------------------------------------

#if __KIND == __SINGLE_PRECISION
         CALL ppm_mg_prolong_2d_vec_s(1,info)
#elif __KIND == __DOUBLE_PRECISION
         CALL ppm_mg_prolong_2d_vec_d(1,info)
#endif   
         !----------------------------------------------------------------------
         !UPDATE THE FUNCTION
         !----------------------------------------------------------------------
         DO isub=1,nsubs
            tuc=>mgfield(isub,mlev)%uc
            DO j=start(2,isub,1),istop(2,isub,1)   
               DO i=start(1,isub,1),istop(1,isub,1)
                DO ilda=1,vecdim
                     u(ilda,i,j,isub)=tuc(ilda,i,j)
                ENDDO 
               ENDDO
            ENDDO
         ENDDO

         !----------------------------------------------------------------------
         !DO the final sweeps
           !--------------------------------------------------------------------
         CALL ppm_mg_smooth_vec(u,f,iterf,1,c1,c2,c3,info) 
         CALL ppm_mg_res_vec(u,f,c1,c2,c3,c4,E,info)   

#ifdef USE_MPI
         CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
         Eu=gEu        
#else
         Eu=E
#endif

#elif __MESH_DIM == __3D



         c1 = 1.0_MK/(2.0_MK*(rdx2+rdy2+rdz2))
         c2 = rdx2
         c3 = rdy2
         c4 = rdz2 
         c5 = 2.0_MK*c2+2.0_MK*c3+2.0_MK*c4


         CALL ppm_mg_smooth_vec(u,f,itera,1,c1,c2,c3,c4,info) 
         !----------------------------------------------------------------------
         ! Compute residual
         !----------------------------------------------------------------------

         CALL ppm_mg_res_vec(u,f,c1,c2,c3,c4,c5,E,info)   

#ifdef USE_MPI
         CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
         E=gEu        
#endif

         IF (l_print) THEN 
          WRITE(cbuf,*) 'Eu:',E
          CALL PPM_WRITE(ppm_rank,'mg_solv',cbuf,info)
         ENDIF

          !---------------------------------------------------------------------
         !Initiation of the function correction. (We start on purpose with lev=2)
         !----------------------------------------------------------------------
         DO mlev=2,maxlev
            DO isub=1,nsubs
               tuc=>mgfield(isub,mlev)%uc
               DO k=start(3,isub,mlev),istop(3,isub,mlev) 
                  DO j=start(2,isub,mlev),istop(2,isub,mlev)
                     DO i=start(1,isub,mlev),istop(1,isub,mlev)
#ifdef __VECTOR
                           tuc(1,i,j,k)=0.0_MK
                           tuc(2,i,j,k)=0.0_MK
                           tuc(3,i,j,k)=0.0_MK
#else
                      DO ilda=1,vecdim 
                           tuc(ilda,i,j,k)=0.0_MK
                      ENDDO
#endif
                     ENDDO
                 ENDDO
               ENDDO
            ENDDO
         ENDDO
         !----------------------------------------------------------------------
         !CALL THE MULTIGRID TO DO NICE STUFF TO OUR FUNCTION
         !----------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
         CALL ppm_mg_core_3d_vec_s(2,iter1,iter2,info)
#elif __KIND == __DOUBLE_PRECISION
         CALL ppm_mg_core_3d_vec_d(2,iter1,iter2,info)
#endif
         !----------------------------------------------------------------------
         !PROLONG the solution to the finest grid
         !----------------------------------------------------------------------

#if __KIND == __SINGLE_PRECISION
         CALL ppm_mg_prolong_3d_vec_s(1,info)
#elif __KIND == __DOUBLE_PRECISION
         CALL ppm_mg_prolong_3d_vec_d(1,info)
#endif

         !----------------------------------------------------------------------
         !UPDATE THE FUNCTION
         !----------------------------------------------------------------------
         DO isub=1,nsubs
            tuc=>mgfield(isub,mlev)%uc
            DO k=start(3,isub,1),istop(3,isub,1)
               DO j=start(2,isub,1),istop(2,isub,1)
                  DO i=start(1,isub,1),istop(1,isub,1)
#ifdef __VECTOR
                        u(1,i,j,k,isub)=tuc(1,i,j,k)
                        u(2,i,j,k,isub)=tuc(2,i,j,k)
                        u(3,i,j,k,isub)=tuc(3,i,j,k)
#else
                   DO ilda=1,vecdim
                        u(ilda,i,j,k,isub)=tuc(ilda,i,j,k)
                   ENDDO
#endif
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         !----------------------------------------------------------------------
         !DO the final sweeps
           !--------------------------------------------------------------------
         CALL ppm_mg_smooth_vec(u,f,iterf,1,c1,c2,c3,c4,info) 
         CALL ppm_mg_res_vec(u,f,c1,c2,c3,c4,c5,E,info)   


#ifdef USE_MPI
         CALL MPI_AllReduce(E,gEu,1,MPI_PREC,MPI_MAX,ppm_comm,info)
         Eu=gEu
#else
         Eu=E
#endif
#endif
#endif

         !----------------------------------------------------------------------
         !  Return 
         !----------------------------------------------------------------------
 9999    CONTINUE
         CALL substop('ppm_mg_solv',t0,info)
         RETURN
#if    __DIM == __SFIELD
#if    __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_solv_2d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_solv_2d_sca_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_solv_3d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_solv_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if    __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_solv_2d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_solv_2d_vec_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_solv_3d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_solv_3d_vec_d
#endif
#endif
#endif
