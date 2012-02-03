#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!-------------------------------------------------------------------------------
!  Subroutine   :            ppm_mg_res 
!-------------------------------------------------------------------------------
!  Purpose      : In this routine we compute the residual in each level
!            
!                  
!  
!  Input       :
!  Input/output :
! 
!  Output       : info        (I) return status. 0 upon success
!
!  Remarks      :
!
!  References   :
!
!  Revisions    :
!-------------------------------------------------------------------------------
!  $Log: ppm_mg_res_fine.f,v $
!  Revision 1.7  2006/07/21 11:30:56  kotsalie
!  FRIDAY
!
!  Revision 1.5  2006/02/08 19:56:02  kotsalie
!  fixed multiple domains
!
!  Revision 1.4  2005/12/08 12:44:45  kotsalie
!  commiting dirichlet
!
!  Revision 1.3  2005/03/14 13:22:58  kotsalie
!  COMMITED THE VECTOR CASE.IT IS FOR LDA=3
!
!  Revision 1.2  2004/09/28 14:07:19  kotsalie
!  Changes concerning 4th order
!
!  Revision 1.1  2004/09/22 18:46:21  kotsalie
!  MG new version
!
  !-----------------------------------------------------------------------------
!  Parallel Particle Mesh Library (PPM)
!  Institute of Computational Science
!  ETH Zentrum, Hirschengraben 84
!  CH-8092 Zurich, Switzerland
 !------------------------------------------------------------------------------

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_res_fine_2D_sca_s(u,f,c1,c2,c3,c4,E,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_res_fine_2D_sca_d(u,f,c1,c2,c3,c4,E,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_res_fine_3D_sca_s(u,f,c1,c2,c3,c4,c5,E,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_res_fine_3D_sca_d(u,f,c1,c2,c3,c4,c5,E,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_res_fine_2D_vec_s(u,f,c1,c2,c3,c4,E,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_res_fine_2D_vec_d(u,f,c1,c2,c3,c4,E,info)
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_res_fine_3D_vec_s(u,f,c1,c2,c3,c4,c5,E,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_res_fine_3D_vec_d(u,f,c1,c2,c3,c4,c5,E,info)
#endif
#endif
#endif

         !----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
#include "ppm_define.h"

            !-------------------------------------------------------------------
        !  Modules 
        !-----------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_data_mg
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_error
        USE ppm_module_alloc
        USE ppm_module_data_mesh



        IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#else
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
            !-------------------------------------------------------------------
        !  Arguments     
        !-----------------------------------------------------------------------
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
#if  __MESH_DIM == __2D
        REAL(MK),                  INTENT(IN)      ::  c1,c2,c3,c4 
#elif __MESH_DIM == __3D
        REAL(MK),                  INTENT(IN)      ::  c1,c2,c3,c4,c5 
#endif
        REAL(MK),                  INTENT(OUT)     ::  E
        INTEGER,                   INTENT(INOUT)   ::  info
          !---------------------------------------------------------------------
        !  Local variables 
        !-----------------------------------------------------------------------
        CHARACTER(LEN=256) :: cbuf
        INTEGER                                    ::  i,j,isub,color
        INTEGER                                    ::  ilda,isweep,count
        INTEGER                                    ::  aa,bb,cc,dd,ee,gg
        REAL(MK)                                   ::  c11,c22,c33,c44,c55 
        INTEGER                                    ::  k,idom
        REAL(MK)                                   ::  x,y
        REAL(MK)                                   ::  res
#if __MESH_DIM == __2D
        INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
        INTEGER,DIMENSION(3)                       ::  ldl3,ldu3
#endif
#if __MESH_DIM == __3D
        INTEGER,DIMENSION(5)                       ::  ldl5,ldu5
        INTEGER,DIMENSION(4)                       ::  ldl4,ldu4
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

        !-----------------------------------------------------------------------
        !Externals
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !Initialize
        !-----------------------------------------------------------------------

        CALL substart('ppm_mg_res',t0,info)
         

        !-----------------------------------------------------------------------
        !  Check arguments
        !-----------------------------------------------------------------------
        IF (ppm_debug .GT. 0) THEN
          IF (c1.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_res',  &
     &            'Factor c1 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (c2.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_res',  &
     &            'Factor c2 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (c3.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_res',  &
     &            'Factor c3 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (c4.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_res',  &
     &            'Factor c4 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#if __MESH_DIM == __3D
          IF (c5.LE.0.0_MK) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_res',  &
     &            'Factor c5 must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
#endif
        ENDIF
        !-----------------------------------------------------------------------
        !Definition of necessary variables and allocation of arrays
        !-----------------------------------------------------------------------
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


#if  __DIM == __SFIELD
#if  __MESH_DIM == __2D

        !-----------------------------------------------------------------------
        !Implementation
         !----------------------------------------------------------------------
          E =-HUGE(E)
          DO isub=1,nsubs
            DO j=start(2,isub,1),istop(2,isub,1)
               DO i=start(1,isub,1),istop(1,isub,1)
                     res  = (u(i-1,j,isub)+u(i+1,j,isub))*c2 + &
     &                      (u(i,j-1,isub)+u(i,j+1,isub))*c3 - &
     &                       u(i,j,isub)*c4-f(i,j,isub)

                     E   = MAX(E,abs(res))
                     mgfield(isub,1)%err(i,j)=-res
                     mgfield(isub,1)%uc(i,j)=u(i,j,isub)
               ENDDO
            ENDDO
          ENDDO

#elif __MESH_DIM == __3D

        !-----------------------------------------------------------------------
        !Implementation
         !----------------------------------------------------------------------

                 E =-HUGE(E)
                 DO isub=1,nsubs
                   aa=0
                   bb=0
                   cc=0
                   dd=0
                   ee=0
                   gg=0

                IF (.NOT.lperiodic) THEN
                 DO iface=1,6
                    IF (bcdef_sca(isub,iface).EQ.ppm_param_bcdef_periodic) THEN
                   !DO NOTHING
                    ELSEIF (bcdef_sca(isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                      IF (iface.EQ.1) THEN
                        aa=1
                      ELSEIF (iface.EQ.2) THEN
                        bb=1
                      ELSEIF (iface.EQ.3) THEN
                       cc=1
                      ELSEIF (iface.EQ.4) THEN
                       dd=1
                      ELSEIF (iface.EQ.5) Then
                       ee=1
                      ELSEIF (iface.EQ.6) Then
                       gg=1
                             ENDIF
                     ENDIF 
                    ENDDO !iface
               ENDIF !periodic

           DO k=start(3,isub,1)+ee,istop(3,isub,1)-gg
              DO j=start(2,isub,1)+cc,istop(2,isub,1)-dd
                DO i=start(1,isub,1)+aa,istop(1,isub,1)-bb
                       res  = (u(i-1,j,k,isub)+u(i+1,j,k,isub))*c2 + &
     &                        (u(i,j-1,k,isub)+u(i,j+1,k,isub))*c3 + &
     &                        (u(i,j,k-1,isub)+u(i,j,k+1,isub))*c4 - &
     &                         u(i,j,k,isub)*c5-f(i,j,k,isub)
                       E   = MAX(E,abs(res))
                       mgfield(isub,1)%err(i,j,k)=-res
                       mgfield(isub,1)%uc(i,j,k)=u(i,j,k,isub)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO



#endif
#elif __DIM == __VFIELD
#if  __MESH_DIM == __2D

        !-----------------------------------------------------------------------
        !Implementation
         !----------------------------------------------------------------------
        E =-HUGE(E)
        DO isub=1,nsubs
           DO j=start(2,isub,1),istop(2,isub,1)
              DO i=start(1,isub,1),istop(1,isub,1)
               DO ilda=1,vecdim
                    res  = (u(ilda,i-1,j,isub)+u(ilda,i+1,j,isub))*c2 + &
     &                     (u(ilda,i,j-1,isub)+u(ilda,i,j+1,isub))*c3 - &
     &                      u(ilda,i,j,isub)*c4-f(ilda,i,j,isub)

                    E   = MAX(E,abs(res))
                    mgfield(isub,1)%err(ilda,i,j)=-res
                    mgfield(isub,1)%uc(ilda,i,j)=u(ilda,i,j,isub)
               ENDDO
              ENDDO
           ENDDO
        ENDDO
        


#elif __MESH_DIM == __3D

        !-----------------------------------------------------------------------
        !Implementation
         !----------------------------------------------------------------------

        

                 E =-HUGE(E)
                 DO isub=1,nsubs
                   aa=0
                   bb=0
                   cc=0
                   dd=0
                   ee=0
                   gg=0
                   DO ilda=1,vecdim

                IF (.NOT.lperiodic) THEN
                 DO iface=1,6
                    IF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_periodic) THEN
                   !DO NOTHING
                    ELSEIF (bcdef_vec(ilda,isub,iface).EQ.ppm_param_bcdef_dirichlet) THEN
                      IF (iface.EQ.1) THEN
                        aa=1
                      ELSEIF (iface.EQ.2) THEN
                        bb=1
                      ELSEIF (iface.EQ.3) THEN
                       cc=1
                      ELSEIF (iface.EQ.4) THEN
                       dd=1
                      ELSEIF (iface.EQ.5) Then
                       ee=1
                      ELSEIF (iface.EQ.6) Then
                       gg=1
                             ENDIF
                     ENDIF 
                    ENDDO !iface
               ENDIF !periodic
              ENDDO

           DO k=start(3,isub,1)+ee,istop(3,isub,1)-gg
              DO j=start(2,isub,1)+cc,istop(2,isub,1)-dd
                DO i=start(1,isub,1)+aa,istop(1,isub,1)-bb
#ifdef __VECTOR
                       res  = (u(1,i-1,j,k,isub)+u(1,i+1,j,k,isub))*c2 +&
     &                        (u(1,i,j-1,k,isub)+u(1,i,j+1,k,isub))*c3 +&
     &                        (u(1,i,j,k-1,isub)+u(1,i,j,k+1,isub))*c4 -&
     &                         u(1,i,j,k,isub)*c5-f(1,i,j,k,isub)

                       E   = MAX(E,abs(res))
                       mgfield(isub,1)%err(1,i,j,k)=-res
                       mgfield(isub,1)%uc(1,i,j,k)=u(1,i,j,k,isub)

                       res  = (u(2,i-1,j,k,isub)+u(2,i+1,j,k,isub))*c2 +&
     &                        (u(2,i,j-1,k,isub)+u(2,i,j+1,k,isub))*c3 +&
     &                        (u(2,i,j,k-1,isub)+u(2,i,j,k+1,isub))*c4 -&
     &                         u(2,i,j,k,isub)*c5-f(2,i,j,k,isub)

                       E   = MAX(E,abs(res))
                       mgfield(isub,1)%err(2,i,j,k)=-res
                       mgfield(isub,1)%uc(2,i,j,k)=u(2,i,j,k,isub)


                       res  = (u(3,i-1,j,k,isub)+u(3,i+1,j,k,isub))*c2 +&
     &                        (u(3,i,j-1,k,isub)+u(3,i,j+1,k,isub))*c3 +&
     &                        (u(3,i,j,k-1,isub)+u(3,i,j,k+1,isub))*c4 -&
     &                         u(3,i,j,k,isub)*c5-f(3,i,j,k,isub)

                       E   = MAX(E,abs(res))
                       mgfield(isub,1)%err(3,i,j,k)=-res
                       mgfield(isub,1)%uc(3,i,j,k)=u(3,i,j,k,isub)

#else
                 DO ilda=1,vecdim
                       res  = (u(ilda,i-1,j,k,isub)+u(ilda,i+1,j,k,isub))*c2 +&
     &                        (u(ilda,i,j-1,k,isub)+u(ilda,i,j+1,k,isub))*c3 +&
     &                        (u(ilda,i,j,k-1,isub)+u(ilda,i,j,k+1,isub))*c4 -&
     &                         u(ilda,i,j,k,isub)*c5-f(ilda,i,j,k,isub)

                       E   = MAX(E,abs(res))
                       mgfield(isub,1)%err(ilda,i,j,k)=-res
                       mgfield(isub,1)%uc(ilda,i,j,k)=u(ilda,i,j,k,isub)
                  ENDDO
#endif
                 ENDDO
              ENDDO
           ENDDO
        ENDDO

#endif
#endif


         !----------------------------------------------------------------------
        !  Return 
        !-----------------------------------------------------------------------
9999    CONTINUE
        CALL substop('ppm_mg_res',t0,info)
        RETURN
#if __DIM == __SFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_res_fine_2D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_res_fine_2D_sca_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_res_fine_3D_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_res_fine_3D_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if   __MESH_DIM   == __2D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_res_fine_2D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_res_fine_2D_vec_d
#endif
#elif __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_res_fine_3D_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_res_fine_3D_vec_d
#endif
#endif
#endif




