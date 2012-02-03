#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


 !------------------------------------------------------------------------------
 !  Subroutine   :            ppm_mg_smooth_coarse    
 !------------------------------------------------------------------------------
 !  Purpose      : In this routine we compute the corrections for
 !                 the function based on the Gauss-Seidel iteration
 !                  
 !  
 !  Input        : nsweep      (I) number of iterations(sweeps)
 !  Input/output :
 ! 
 !  Output       : info        (I) return status. 0 upon success
 !
 !  Remarks      :
 !
 !  References   :
 !
 !  Revisions    :
 !------------------------------------------------------------------------------
 !  $Log: ppm_mg_smooth_coarse.f,v $
 !  Revision 1.14  2006/07/21 11:30:55  kotsalie
 !  FRIDAY
 !
 !  Revision 1.12  2006/02/08 19:55:05  kotsalie
 !  fixed multiple subdomains
 !
 !  Revision 1.11  2006/02/02 17:59:45  michaebe
 !  corrected a bug in the log comment
 !
 !  Revision 1.10  2006/02/02 16:33:19  kotsalie
 !  corrected for mixed bc''s
 !
 !  Revision 1.9  2005/12/08 12:44:46  kotsalie
 !  commiting dirichlet
 !
 !  Revision 1.8  2005/03/14 13:27:32  kotsalie
 !  COMMITED THE VECTOR CASE. IT IS FOR LDA=3
 !
 !  Revision 1.7  2005/01/04 09:45:29  kotsalie
 !  ghostsize=2
 !
 !  Revision 1.6  2004/11/05 18:09:49  kotsalie
 !  FINAL FEATURE BEFORE TEST.I DO NOT USE MASKS
 !
 !  Revision 1.4  2004/10/29 15:59:31  kotsalie
 !  RED BLACK SOR
 !
 !  Revision 1.3  2004/09/28 14:05:31  kotsalie
 !  Changes concernig 4th order finite differences
 !
 !  Revision 1.2  2004/09/23 12:16:49  kotsalie
 !  Added USE statement
 !
 !  Revision 1.1  2004/09/22 18:42:39  kotsalie
 !  MG new version
 !
 !
   !----------------------------------------------------------------------------
 !  Parallel Particle Mesh Library (PPM)
 !  Institute of Computational Science
 !  ETH Zentrum, Hirschengraben 84
 !  CH-8092 Zurich, Switzerland
  !-----------------------------------------------------------------------------

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_smooth_coarse_2D_sca_s(nsweep,mlev,c1,c2,c3,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_smooth_coarse_2D_sca_d(nsweep,mlev,c1,c2,c3,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_smooth_coarse_3D_sca_s(nsweep,mlev,c1,c2,c3,c4,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_smooth_coarse_3D_sca_d(nsweep,mlev,c1,c2,c3,c4,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_smooth_coarse_2D_vec_s(nsweep,mlev,c1,c2,c3,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_smooth_coarse_2D_vec_d(nsweep,mlev,c1,c2,c3,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_smooth_coarse_3D_vec_s(nsweep,mlev,c1,c2,c3,c4,info)
#elif  __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_smooth_coarse_3D_vec_d(nsweep,mlev,c1,c2,c3,c4,info)
#endif
#endif
#endif

          !---------------------------------------------------------------------
         !  Includes
         !----------------------------------------------------------------------
#include "ppm_define.h"

             !------------------------------------------------------------------
         !  Modules 
         !----------------------------------------------------------------------
         USE ppm_module_data
         USE ppm_module_data_mg
         USE ppm_module_substart
         USE ppm_module_substop
         USE ppm_module_error
         USE ppm_module_alloc
         USE ppm_module_map_field_ghost
         USE ppm_module_data_mesh
         USE ppm_module_write



         IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
         INTEGER, PARAMETER :: MK = ppm_kind_single
#else
         INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
         !------------------------------------------------------------------
         !  Arguments     
         !----------------------------------------------------------------------
         INTEGER,                   INTENT(IN)      ::  nsweep
         INTEGER,                   INTENT(IN)      ::  mlev
#if  __MESH_DIM == __2D
         REAL(MK),                  INTENT(IN)      ::  c1,c2,c3 
#elif __MESH_DIM == __3D
         REAL(MK),                  INTENT(IN)      ::  c1,c2,c3,c4 
#endif
         INTEGER,                   INTENT(INOUT)   ::  info
           !--------------------------------------------------------------------
         !  Local variables 
         !----------------------------------------------------------------------
         CHARACTER(LEN=256) :: cbuf
         INTEGER                                    ::  i,j,isub,color
         INTEGER,DIMENSION(:),POINTER               ::  a,b,c,d,e,g 
         REAL(MK)                                   ::  c11,c22,c33,c44 
         INTEGER                                    ::  ilda,isweep,count
         INTEGER                                    ::  k,idom
         REAL(MK)                                   ::  x,y,dx,dy
         REAL(MK)                                   ::  omega
         INTEGER,DIMENSION(1)                       ::  ldu1,ldl1
#if __MESH_DIM == __2D
         INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
         INTEGER,DIMENSION(3)                       ::  ldl3,ldu3
#endif
#if __MESH_DIM == __3D
         INTEGER,DIMENSION(5)                       ::  ldl5,ldu5
         INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
         REAL(MK)                                   ::  dz
#endif
         INTEGER                                    ::  iopt,iface,topoid
         REAL(MK)                                   ::  t0
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
         REAL(MK),DIMENSION(:,:,:),POINTER :: oldu
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:),POINTER :: oldu
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
         REAL(MK),DIMENSION(:,:,:,:),POINTER :: oldu  
#elif __MESH_DIM == __3D
         REAL(MK),DIMENSION(:,:,:,:,:),POINTER :: oldu
#endif
#endif

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
      REAL(MK) :: moldu
#elif __MESH_DIM == __3D
      REAL(MK) :: moldu
#endif
#elif  __DIM == __VFIELD
#if __MESH_DIM == __2D
      REAL(MK),DIMENSION(:),POINTER :: moldu
#elif __MESH_DIM == __3D
      REAL(MK),DIMENSION(:),POINTER :: moldu
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


#if __KIND == __SINGLE_PRECISION
       omega=omega_s
       dx=dx_s
       dy=dy_s
#if __MESH_DIM == __3D
       dz=dz_s
#endif
#elif __KIND == __DOUBLE_PRECISION
       omega=omega_d
       dx=dx_d
       dy=dy_d
#if __MESH_DIM == __3D
       dz=dz_d
#endif
#endif

         !----------------------------------------------------------------------
         !Externals
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         !Initialize
         !----------------------------------------------------------------------

         CALL substart('ppm_mg_smooth_coarse',t0,info)
         IF (l_print) THEN 
          WRITE (cbuf,*) 'SMOOTHER entering ','mlev:',mlev
          CALL PPM_WRITE(ppm_rank,'mg_smooth',cbuf,info)
         ENDIF

         !----------------------------------------------------------------------
         !  Check arguments
         !----------------------------------------------------------------------
         IF (ppm_debug .GT. 0) THEN
           IF (nsweep.LT.1) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'nsweep must be >=1',__LINE__,info)
               GOTO 9999
           ENDIF
           IF (mlev.LE.1) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'level must be >1',__LINE__,info)
               GOTO 9999
           ENDIF
           IF (c1.LE.0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'Factor c1 must be >0',__LINE__,info)
               GOTO 9999
           ENDIF
           IF (c2.LE.0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'Factor c2 must be >0',__LINE__,info)
               GOTO 9999
           ENDIF
           IF (c3.LE.0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'Factor c3 must be >0',__LINE__,info)
               GOTO 9999
           ENDIF
#if __MESH_DIM == __3D
           IF (c4.LE.0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_mg_smooth_coarse',  &
      &            'Factor c4 must be >0',__LINE__,info)
               GOTO 9999
          ENDIF
#endif
         ENDIF
         !----------------------------------------------------------------------
         !Definition of necessary variables and allocation of arrays
         !----------------------------------------------------------------------
         topoid=ppm_field_topoid


#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         mgfield=>mgfield_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
         mgfield=>mgfield_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
         mgfield=>mgfield_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
         mgfield=>mgfield_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
         mgfield=>mgfield_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
         mgfield=>mgfield_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
         mgfield=>mgfield_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
         mgfield=>mgfield_3d_vec_d
#endif
#endif
#endif

             iopt = ppm_param_alloc_fit
             ldl1(1) = 1
             ldu1(1) = nsubs
             CALL ppm_alloc(a,ldl1,ldu1,iopt,info)
             CALL ppm_alloc(b,ldl1,ldu1,iopt,info)
             CALL ppm_alloc(c,ldl1,ldu1,iopt,info)
             CALL ppm_alloc(d,ldl1,ldu1,iopt,info)
             CALL ppm_alloc(e,ldl1,ldu1,iopt,info)
             CALL ppm_alloc(g,ldl1,ldu1,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'GSsolv',    &
       &                       'a',__LINE__,info)
             GOTO 9999
             ENDIF

#if  __DIM == __SFIELD
#if  __MESH_DIM == __2D

         !----------------------------------------------------------------------
         !Implementation
          !---------------------------------------------------------------------

             iopt = ppm_param_alloc_fit
             ldl3(1) = 1-ghostsize(1)
             ldl3(2) = 1-ghostsize(2)
             ldl3(3) = 1
             ldu3(1) = max_node(1,mlev)+ghostsize(1)
             ldu3(2) = max_node(2,mlev)+ghostsize(2)
             ldu3(3) = nsubs
             CALL ppm_alloc(uc_dummy,ldl3,ldu3,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'GSsolv',    &
       &                       'uc_dummy',__LINE__,info)
             GOTO 9999
             ENDIF


         DO isweep=1,nsweep
            DO color=0,1

               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc
                  uc_dummy(:,:,isub)=tuc(:,:)
               ENDDO!DO isub 
               !----------------------------------------------------------------
               !Communicate
               !----------------------------------------------------------------


               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                    ghostsize,ppm_param_map_ghost_get,info) 
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_push,info) 
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_send,info) 
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                          ghostsize,ppm_param_map_pop,info) 



               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc
                           tuc(:,:)=uc_dummy(&
      &                         :,:,isub)
                DO j=start(2,isub,mlev),istop(2,isub,mlev)
                   DO i=start(1,isub,mlev)+mod(j+color,2),&
                         &istop(1,isub,mlev),2
                           tuc(i,j) = c1*(&
      &                                   (tuc(i-1,j)+ &
      &                                tuc(i+1,j))*c2 + &
      &                                 (tuc(i,j-1)+&
      &                                  tuc(i,j+1))*c3-&
      &                                         mgfield(isub,mlev)%fc(i,j))
                     ENDDO
                  ENDDO
               ENDDO!isub

               IF (isweep.EQ.nsweep) THEN   
                IF (color.EQ.1) THEN
                 DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc
                  uc_dummy(:,:,isub)=tuc(:,:) 
                 ENDDO   
                ENDIF
               ENDIF
              ENDDO!DO color   

              IF (isweep.EQ.nsweep) THEN

               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                    ghostsize,ppm_param_map_ghost_get,info)
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_push,info)
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_send,info)
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                          ghostsize,ppm_param_map_pop,info)

               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc
                           tuc(:,:)=uc_dummy(&
      &                         :,:,isub)
               ENDDO  
             ENDIF


            ENDDO!DO nsweep



             iopt = ppm_param_dealloc
             ldl3(1) = 1-ghostsize(1)
             ldl3(2) = 1-ghostsize(2)
             ldl3(3) = 1
             ldu3(1) = max_node(1,mlev)+ghostsize(1)
             ldu3(2) = max_node(2,mlev)+ghostsize(2)
             ldu3(3) = nsubs
             CALL ppm_alloc(uc_dummy,ldl3,ldu3,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'GSsolv',    &
       &                       'uc_dummy',__LINE__,info)
             GOTO 9999
             ENDIF

#elif __MESH_DIM == __3D

         !----------------------------------------------------------------------
         !Implementation
          !---------------------------------------------------------------------

             iopt = ppm_param_alloc_fit
             ldl4(1) = 1-ghostsize(1)
             ldl4(2) = 1-ghostsize(2)
             ldl4(3) = 1-ghostsize(3)
             ldl4(4) = 1
             ldu4(1) = max_node(1,mlev)+ghostsize(1)
             ldu4(2) = max_node(2,mlev)+ghostsize(2)
             ldu4(3) = max_node(3,mlev)+ghostsize(3)
             ldu4(4) = nsubs
             CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'GSsolv',    &
       &                       'uc_dummy',__LINE__,info)
             GOTO 9999
             ENDIF


         DO isweep=1,nsweep 
            DO color=0,1
               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc  
                    DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                     DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                      DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                        uc_dummy(i,j,k,isub)=tuc(i,j,k)
                      ENDDO
                     ENDDO
                    ENDDO
               ENDDO!DO isub 


               !----------------------------------------------------------------
               !Communicate
               !----------------------------------------------------------------


               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                    ghostsize,ppm_param_map_ghost_get,info) 
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_push,info) 
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_send,info) 
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                          ghostsize,ppm_param_map_pop,info) 

                 a=0
                 b=0
                 c=0
                 d=0
                 e=0
                 g=0

               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc  
                    DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                     DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                      DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                          tuc(i,j,k)=uc_dummy(i,j,k,isub)
                      ENDDO
                     ENDDO
                   ENDDO
                 IF (.NOT.lperiodic) THEN
                  DO iface=1,6
                   IF (bcdef_sca(isub,iface).EQ.ppm_param_bcdef_periodic) THEN
                    !DO NOTHING 
                   ELSEIF (bcdef_sca(isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN

                     IF (iface.EQ.1) THEN
                        a(isub)=1
                        IF (bcdef_sca(isub,2).EQ.0) THEN
                         b(isub)=-1  
                        ENDIF 
                       i=1
                        DO j=1,max_node(2,mlev)
                         DO k=1,max_node(3,mlev)
                           tuc(i,j,k)=0.0_MK
                         enddo
                        ENDDO
                     ELSEIF (iface.EQ.2) THEN
                        b(isub)=1
                        IF (bcdef_sca(isub,1).EQ.0) THEN
                         a(isub)=-1  
                        ENDIF 
                       i=max_node(1,mlev)
                        DO j=1,max_node(2,mlev)
                         DO k=1,max_node(3,mlev)

                          tuc(i,j,k)=0.0_MK
                         ENDDO
                        enddo
                     ELSEIF (iface.EQ.3) THEN
                       c(isub)=1
                        IF (bcdef_sca(isub,4).EQ.0) THEN
                         d(isub)=-1  
                        ENDIF 
                       j=1
                        DO i=1,max_node(1,mlev)
                         Do k=1,max_node(3,mlev)
                          tuc(i,j,k)=0.0_MK
                         enddo
                        ENDDO
                     ELSEIF (iface.EQ.4) THEN
                       d(isub)=1
                        IF (bcdef_sca(isub,3).EQ.0) THEN
                         c(isub)=-1  
                        ENDIF 
                       j=max_node(2,mlev)
                        DO i=1,max_node(1,mlev)
                         Do k=1,max_node(3,mlev)
                          tuc(i,j,k)=0.0_MK
                         enddo
                        ENDDO
                     ELSEIF (iface.EQ.5) Then
                       e(isub)=1
                        IF (bcdef_sca(isub,6).EQ.0) THEN
                         g(isub)=-1  
                        ENDIF 
                       k=1
                        DO i=1,max_node(1,mlev)
                         Do j=1,max_node(2,mlev)
                          tuc(i,j,k)=0.0_MK
                         enddo
                        ENDDO
                              ELSEIF (iface.EQ.6) Then
                       g(isub)=1
                        IF (bcdef_sca(isub,5).EQ.0) THEN
                         e(isub)=-1  
                        ENDIF 
                        k=max_node(3,mlev)
                                DO i=1,max_node(1,mlev) 
                                 Do j=1,max_node(2,mlev)
                          tuc(i,j,k)=0.0_MK
                                     enddo
                        ENDDO
                               endif                  
                  ENDIF 
                 ENDDO!iface 
                ENDIF 
                  DO k=start(3,isub,mlev)+e(isub),istop(3,isub,mlev)-g(isub) 
                     DO j=start(2,isub,mlev)+c(isub),istop(2,isub,mlev)-d(isub)
                        DO i=start(1,isub,mlev)+mod(j+k+color,2)+a(isub),&
                            & istop(1,isub,mlev)-b(isub)-mod(j+k+color,2),2

                         IF ((i.GE.1.AND.i.LE.max_node(1,mlev)).AND. &
                              (j.GE.1.AND.j.LE.max_node(2,mlev)).AND.(k.GE.1.AND.k.LE.max_node(3,mlev))) THEN
                             moldu=tuc(i,j,k)

                              tuc(i,j,k) = moldu+&
      &                             omega*(&
      &                             c1*((tuc(i-1,j,k)+ &
      &                            tuc(i+1,j,k))*c2 + &
      &                                 (tuc(i,j-1,k)+&
      &                            tuc(i,j+1,k))*c3 + &
      &                           (tuc(i,j,k-1)+&
      &                            tuc(i,j,k+1))*c4 - &
      &                                    mgfield(isub,mlev)%fc(i,j,k))&
      &                            -moldu)
                         ENDIF
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO!isubs   

                   IF (isweep.EQ.nsweep) THEN  
                     IF (color.EQ.1) THEN
                      DO isub=1,nsubs
                       tuc=>mgfield(isub,mlev)%uc  
                       DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                        DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                         DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                          uc_dummy(i,j,k,isub)=tuc(i,j,k)
                         ENDDO
                        ENDDO
                       ENDDO

                     ENDDO!isub
                    ENDIF
                  ENDIF

           ENDDO!DO color

               IF (isweep.EQ.nsweep) THEN


               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                    ghostsize,ppm_param_map_ghost_get,info) 
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_push,info) 
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_send,info) 
               CALL ppm_map_field_ghost(uc_dummy,topoid,mesh_id_g(mlev),&
      &                          ghostsize,ppm_param_map_pop,info) 

               ENDIF

               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc  

                    DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                     DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                      DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                          tuc(i,j,k)=uc_dummy(i,j,k,isub)
                      ENDDO
                     ENDDO
                   ENDDO

               ENDDO!isub
         ENDDO!Do isweep

             iopt = ppm_param_dealloc
             ldl4(1) = 1-ghostsize(1)
             ldl4(2) = 1-ghostsize(2)
             ldl4(3) = 1-ghostsize(3)
             ldl4(4) = 1
             ldu4(1) = max_node(1,mlev)+ghostsize(1)
             ldu4(2) = max_node(2,mlev)+ghostsize(2)
             ldu4(3) = max_node(3,mlev)+ghostsize(3)
             ldu4(4) = nsubs
             CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'GSsolv',    &
       &                       'uc_dummy',__LINE__,info)
             GOTO 9999
             ENDIF
#endif
#elif __DIM == __VFIELD
#if  __MESH_DIM == __2D

         !----------------------------------------------------------------------
         !Implementation
          !---------------------------------------------------------------------

             iopt = ppm_param_alloc_fit
             ldl4(1) = 1
             ldl4(2) = 1-ghostsize(1)
             ldl4(3) = 1-ghostsize(2)
             ldl4(4) = 1
             ldu4(1) = vecdim
             ldu4(2) = max_node(1,mlev)+ghostsize(1)
             ldu4(3) = max_node(2,mlev)+ghostsize(2)
             ldu4(4) = nsubs
             CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'GSsolv',    &
       &                       'uc_dummy',__LINE__,info)
             GOTO 9999
             ENDIF



         DO isweep=1,nsweep
            DO color=0,1
               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc
                  uc_dummy(:,:,:,isub)=tuc(:,:,:)
               ENDDO!DO isub 

               !----------------------------------------------------------------
               !Communicate 
               !----------------------------------------------------------------

               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                    ghostsize,ppm_param_map_ghost_get,info) 
               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_push,info) 
               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_send,info) 
               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                          ghostsize,ppm_param_map_pop,info) 



               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc
                  tuc(:,:,:)=uc_dummy(&
      &                         :,:,:,isub)
                  DO j=start(2,isub,mlev),istop(2,isub,mlev)
                     DO i=start(1,isub,mlev)+mod(j+color,2),istop(1,isub,mlev),2
                      DO ilda=1,vecdim
                           tuc(ilda,i,j) = c1*(&
      &                                   (tuc(ilda,i-1,j)+ &
      &                                tuc(ilda,i+1,j))*c2 + &
      &                                 (tuc(ilda,i,j-1)+&
      &                                  tuc(ilda,i,j+1))*c3-&
      &                                         mgfield(isub,mlev)%fc(ilda,i,j))

                      ENDDO  
                     ENDDO
                  ENDDO
               ENDDO
                    IF (isweep.EQ.nsweep) THEN
                     IF (color.EQ.1) THEN
                      DO isub=1,nsubs
                       tuc=>mgfield(isub,mlev)%uc
                       uc_dummy(:,:,:,isub)=tuc(:,:,:)
                      ENDDO
                     ENDIF
                    ENDIF




            ENDDO!DO color   

              IF (isweep.EQ.nsweep) THEN
               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                    ghostsize,ppm_param_map_ghost_get,info) 
               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_push,info) 
               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_send,info) 
               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                          ghostsize,ppm_param_map_pop,info) 


               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc
                  tuc(:,:,:)=uc_dummy(&
      &                         :,:,:,isub)
               ENDDO
              ENDIF 

         ENDDO!DO nsweep



             iopt = ppm_param_dealloc
             ldl4(1) = 1
             ldl4(2) = 1-ghostsize(1)
             ldl4(3) = 1-ghostsize(2)
             ldl4(4) = 1
             ldu4(1) = vecdim
             ldu4(2) = max_node(1,mlev)+ghostsize(1)
             ldu4(3) = max_node(2,mlev)+ghostsize(2)
             ldu4(4) = nsubs
             CALL ppm_alloc(uc_dummy,ldl4,ldu4,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'GSsolv',    &
       &                       'uc_dummy',__LINE__,info)
             GOTO 9999
             ENDIF
#elif __MESH_DIM == __3D

         !----------------------------------------------------------------------
         !Implementation
          !---------------------------------------------------------------------

             iopt = ppm_param_alloc_fit
             ldl5(1) = 1
             ldl5(2) = 1-ghostsize(1)
             ldl5(3) = 1-ghostsize(2)
             ldl5(4) = 1-ghostsize(3)
             ldl5(5) = 1
             ldu5(1) = vecdim
             ldu5(2) = max_node(1,mlev)+ghostsize(1)
             ldu5(3) = max_node(2,mlev)+ghostsize(2)
             ldu5(4) = max_node(3,mlev)+ghostsize(3)
             ldu5(5) = nsubs
             CALL ppm_alloc(uc_dummy,ldl5,ldu5,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'GSsolv',    &
       &                       'uc_dummy',__LINE__,info)
             GOTO 9999
             ENDIF




            iopt = ppm_param_alloc_fit
            ldu1(1)=vecdim
            CALL ppm_alloc(moldu,ldu1,iopt,info)
            IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'GSsolv',    &
       &                       'moldu',__LINE__,info)
             GOTO 9999
            ENDIF



         DO isweep=1,nsweep 

            DO color=0,1

             DO isub=1,nsubs
                  !-------------------------------------------------------------
                  !Impose boundaries 
                  !-------------------------------------------------------------

                  tuc=>mgfield(isub,mlev)%uc


                    DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                     DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                      DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
#ifdef __VECTOR
                        uc_dummy(1,i,j,k,isub)=tuc(1,i,j,k)
                        uc_dummy(2,i,j,k,isub)=tuc(2,i,j,k)
                        uc_dummy(3,i,j,k,isub)=tuc(3,i,j,k)
#else
                       DO ilda=1,vecdim 
                        uc_dummy(ilda,i,j,k,isub)=tuc(ilda,i,j,k)
                       ENDDO 
#endif
                      ENDDO
                     ENDDO
                    ENDDO 

               ENDDO!DO isub 




               !----------------------------------------------------------------
               !Communicate 
               !----------------------------------------------------------------


               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                    ghostsize,ppm_param_map_ghost_get,info) 
               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_push,info) 
               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                         ghostsize,ppm_param_map_send,info) 
               CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
      &                          ghostsize,ppm_param_map_pop,info) 

                 a=0
                 b=0
                 c=0
                 d=0
                 e=0
                 g=0
               DO isub=1,nsubs
                  tuc=>mgfield(isub,mlev)%uc

                    DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                     DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                      DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
#ifdef __VECTOR

                          tuc(1,i,j,k)=uc_dummy(1,i,j,k,isub)
                          tuc(2,i,j,k)=uc_dummy(2,i,j,k,isub)
                          tuc(3,i,j,k)=uc_dummy(3,i,j,k,isub)

#else

                       DO ilda=1,vecdim 
                          tuc(ilda,i,j,k)=uc_dummy(ilda,i,j,k,isub)
                       ENDDO
#endif
                      ENDDO
                     ENDDO
                   ENDDO

                 Do  ilda=1,vecdim

                  IF (.NOT.lperiodic) THEN

                   DO iface=1,6
                    IF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_periodic) THEN
                     !DO NOTHING
                    ELSEIF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                     IF (iface.EQ.1) THEN
                        a(isub)=1

                        IF (bcdef_vec(ilda,isub,2).EQ.0) THEN
                         b(isub)=-1
                        ENDIF

                       i=1
                        DO j=1,max_node(2,mlev)
                         DO k=1,max_node(3,mlev)
                             tuc(ilda,i,j,k)=0.0_MK
                         enddo
                        ENDDO
                     ELSEIF (iface.EQ.2) THEN
                        b(isub)=1
                        IF (bcdef_vec(ilda,isub,1).EQ.0) THEN
                         a(isub)=-1
                        ENDIF
                       i=max_node(1,mlev)
                        DO j=1,max_node(2,mlev)
                         DO k=1,max_node(3,mlev)
                              tuc(ilda,i,j,k)=0.0_MK
                         ENDDO
                        enddo
                     ELSEIF (iface.EQ.3) THEN
                       c(isub)=1
                        IF (bcdef_vec(ilda,isub,4).EQ.0) THEN
                         d(isub)=-1
                        ENDIF
                       j=1
                        DO i=1,max_node(1,mlev)
                         Do k=1,max_node(3,mlev)
                              tuc(ilda,i,j,k)=0.0_MK

                         enddo
                        ENDDO
                     ELSEIF (iface.EQ.4) THEN
                       d(isub)=1
                        IF (bcdef_vec(ilda,isub,3).EQ.0) THEN
                         c(isub)=-1
                        ENDIF
                       j=max_node(2,mlev)
                        DO i=1,max_node(1,mlev)
                         Do k=1,max_node(3,mlev)
                              tuc(ilda,i,j,k)=0.0_MK
                         enddo
                        ENDDO
                     ELSEIF (iface.EQ.5) Then
                       e(isub)=1
                        IF (bcdef_vec(ilda,isub,6).EQ.0) THEN
                         g(isub)=-1
                        ENDIF
                       k=1
                        DO i=1,max_node(1,mlev)
                         Do j=1,max_node(2,mlev)
                              tuc(ilda,i,j,k)=0.0_MK
                         enddo
                        ENDDO

                      elseif (iface.EQ.6) THEN
                        g(isub)=1
                        IF (bcdef_vec(ilda,isub,5).EQ.0) THEN
                         e(isub)=-1
                        ENDIF
                        k=max_node(3,mlev)
                        DO i=1,max_node(1,mlev)
                         Do j=1,max_node(2,mlev)
                             tuc(ilda,i,j,k)=0.0_MK
                         enddo
                        ENDDO
                      endif
                  ENDIF
                 ENDDO!face
                ENDIF
                ENDDO!ilda
                  DO k=start(3,isub,mlev)+e(isub),istop(3,isub,mlev)-g(isub)  
                     DO j=start(2,isub,mlev)+c(isub),istop(2,isub,mlev)-d(isub)
                        DO i=start(1,isub,mlev)+mod(j+k+color,2)+a(isub),istop(1,isub,mlev)-b(isub)-mod(j+k+color,2),2
                         IF ((i.GE.1.AND.i.LE.max_node(1,mlev)).AND.(j.GE.1.AND.j.LE.max_node(2,mlev)).AND.&
                              (k.GE.1.AND.k.LE.max_node(3,mlev))) THEN
#ifdef __VECTOR

                         moldu(1) = tuc(1,i,j,k)
                         moldu(2) = tuc(2,i,j,k)
                         moldu(3) = tuc(3,i,j,k)
#else
                      do ilda=1,vecdim
                         moldu(ilda) = tuc(ilda,i,j,k)
                      end do
#endif

#ifdef __VECTOR

                              tuc(1,i,j,k) = moldu(1)+&
      &                             omega*(& 
      &                             c1*((tuc(1,i-1,j,k)+ &
      &                            tuc(1,i+1,j,k))*c2 + &
      &                                 (tuc(1,i,j-1,k)+&
      &                            tuc(1,i,j+1,k))*c3 + &
      &                           (tuc(1,i,j,k-1)+&
      &                            tuc(1,i,j,k+1))*c4 - &
      &                            mgfield(isub,mlev)%fc(1,i,j,k))&
      &                            -moldu(1))


                              tuc(2,i,j,k) = moldu(2)+&
      &                             omega*(& 
      &                             c1*((tuc(2,i-1,j,k)+ &
      &                            tuc(2,i+1,j,k))*c2 + &
      &                                 (tuc(2,i,j-1,k)+&
      &                            tuc(2,i,j+1,k))*c3 + &
      &                           (tuc(2,i,j,k-1)+&
      &                            tuc(2,i,j,k+1))*c4 - &
      &                            mgfield(isub,mlev)%fc(2,i,j,k))&
      &                            -moldu(2))

                              tuc(3,i,j,k) = moldu(3)+&
      &                             omega*(& 
      &                             c1*((tuc(3,i-1,j,k)+ &
      &                            tuc(3,i+1,j,k))*c2 + &
      &                                 (tuc(3,i,j-1,k)+&
      &                            tuc(3,i,j+1,k))*c3 + &
      &                           (tuc(3,i,j,k-1)+&
      &                            tuc(3,i,j,k+1))*c4 - &
      &                            mgfield(isub,mlev)%fc(3,i,j,k))&
      &                            -moldu(3))
#else
                      DO ilda=1,vecdim


                              tuc(ilda,i,j,k) = moldu(ilda)+&
      &                             omega*(& 
      &                             c1*((tuc(ilda,i-1,j,k)+ &
      &                            tuc(ilda,i+1,j,k))*c2 + &
      &                                 (tuc(ilda,i,j-1,k)+&
      &                            tuc(ilda,i,j+1,k))*c3 + &
      &                           (tuc(ilda,i,j,k-1)+&
      &                            tuc(ilda,i,j,k+1))*c4 - &
      &                            mgfield(isub,mlev)%fc(ilda,i,j,k))&
      &                            -moldu(ilda))



                         ENDDO!ilda
#endif
                        ENDIF
                        ENDDO!i
                     ENDDO!j
                  ENDDO!k

               ENDDO!isubs   

                   IF (isweep.EQ.nsweep) THEN
                    IF (color.EQ.1) THEN
                     DO isub=1,nsubs

                       tuc=>mgfield(isub,mlev)%uc

                       DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                         DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                           DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                             DO ilda=1,vecdim 
                                 uc_dummy(ilda,i,j,k,isub)=tuc(ilda,i,j,k)
                             ENDDO 
                           ENDDO
                         ENDDO
                       ENDDO 
                     ENDDO!isub   

                    ENDIF
                   ENDIF 


           ENDDO!DO color


          IF (isweep.EQ.nsweep) THEN

           CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
        &             ghostsize,ppm_param_map_ghost_get,info) 
           CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
        &                   ghostsize,ppm_param_map_push,info) 
           CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
        &                   ghostsize,ppm_param_map_send,info) 
           CALL ppm_map_field_ghost(uc_dummy,vecdim,topoid,mesh_id_g(mlev),&
        &                ghostsize,ppm_param_map_pop,info) 


                   DO isub=1,nsubs 
                    tuc=>mgfield(isub,mlev)%uc

                    DO k=1-ghostsize(3),max_node(3,mlev)+ghostsize(3)
                     DO j=1-ghostsize(2),max_node(2,mlev)+ghostsize(2)
                      DO i=1-ghostsize(1),max_node(1,mlev)+ghostsize(1)
                       DO ilda=1,vecdim 
                          tuc(ilda,i,j,k)=uc_dummy(ilda,i,j,k,isub)
                       ENDDO
                      ENDDO
                     ENDDO
                    ENDDO
                   ENDDO
          ENDIF

         ENDDO!Do isweep

             iopt = ppm_param_dealloc
             ldl5(1) = 1
             ldl5(2) = 1-ghostsize(1)
             ldl5(3) = 1-ghostsize(2)
             ldl5(4) = 1-ghostsize(3)
             ldl5(5) = 1
             ldu5(1) = vecdim
             ldu5(2) = max_node(1,mlev)+ghostsize(1)
             ldu5(4) = max_node(2,mlev)+ghostsize(2)
             ldu5(4) = max_node(3,mlev)+ghostsize(3)
             ldu5(5) = nsubs
             CALL ppm_alloc(uc_dummy,ldl5,ldu5,iopt,info)
             IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'GSsolv',    &
       &                       'uc_dummy',__LINE__,info)
             GOTO 9999
             ENDIF

            iopt = ppm_param_dealloc
            ldu1(1)=vecdim
            CALL ppm_alloc(moldu,ldu1,iopt,info)
            IF (info .NE. 0) THEN
             info = ppm_error_fatal
             CALL ppm_error(ppm_err_alloc,'GSsolv',    &
       &                       'moldu',__LINE__,info)
             GOTO 9999
            ENDIF

#endif
#endif


          !---------------------------------------------------------------------
         !  Return 
         !----------------------------------------------------------------------
 9999    CONTINUE
         CALL substop('ppm_mg_smooth_coarse',t0,info)
         RETURN
#if __DIM == __SFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_2D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_2D_sca_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_3D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_3D_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_2D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_2D_vec_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_3D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_smooth_coarse_3D_vec_d
#endif
#endif
#endif




