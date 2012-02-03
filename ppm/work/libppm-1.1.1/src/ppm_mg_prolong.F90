#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

       !------------------------------------------------------------------------
       !  Subroutine   :            ppm_mg_prolong  
       !------------------------------------------------------------------------
       !  Purpose      : In this routine we prolong the corrections from
       !                 coarser to finer levels 
       !                  
       !  
       !  Input        : mlev       (I) current level in V-cycle
       ! 
       !  Input/output :
       ! 
       !  Output       : info       (I) return status. 0 upon success
       !
       !  Remarks      :
       !
       !  References   :
       !
       !  Revisions    :
       !------------------------------------------------------------------------
       !  $Log: ppm_mg_prolong.f,v $
       !  Revision 1.6  2006/07/21 11:30:56  kotsalie
       !  FRIDAY
       !
       !  Revision 1.5  2005/04/21 04:55:03  ivos
       !  fix: corrected misplaced cpp endif statements.
       !
       !  Revision 1.4  2005/03/14 13:17:18  kotsalie
       !  COMMITED THE VECTOR CASE. IT IS FOR LDA=3
       !
       !  Revision 1.3  2004/09/28 14:08:14  kotsalie
       !  *** empty log message ***
       !
       !  Revision 1.2  2004/09/23 12:16:50  kotsalie
       !  Added USE statement
       !
       !  Revision 1.1  2004/09/22 18:36:05  kotsalie
       !  MG new version
       !
        !-----------------------------------------------------------------------
       !  Parallel Particle Mesh Library (PPM)
       !  Institute of Computational Science
       !  ETH Zentrum, Hirschengraben 84
       !  CH-8092 Zurich, Switzerland
        !-----------------------------------------------------------------------


#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_prolong_2d_sca_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_prolong_2d_sca_d(mlev,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_prolong_3d_sca_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_prolong_3d_sca_d(mlev,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_prolong_2d_vec_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_prolong_2d_vec_d(mlev,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
       SUBROUTINE ppm_mg_prolong_3d_vec_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
       SUBROUTINE ppm_mg_prolong_3d_vec_d(mlev,info)
#endif
#endif
#endif

             !------------------------------------------------------------------
         !  Includes
         !----------------------------------------------------------------------
#include "ppm_define.h"

         !----------------------------------------------------------------------
         !  Modules 
         !----------------------------------------------------------------------
         USE ppm_module_data
         USE ppm_module_write
         USE ppm_module_substart         
         USE ppm_module_substop
         USE ppm_module_data_mg
         USE ppm_module_error
         USE ppm_module_alloc
         USE ppm_module_map_field_ghost         

         IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION
         INTEGER, PARAMETER :: MK = ppm_kind_single
#else
         INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
         !----------------------------------------------------------------------
         !  Arguments     
         !----------------------------------------------------------------------
         INTEGER,                   INTENT(IN)      ::  mlev
         INTEGER,                   INTENT(INOUT)   ::  info
              !-----------------------------------------------------------------
         !  Local variables 
         !----------------------------------------------------------------------
         CHARACTER(LEN=256)                         :: cbuf
         INTEGER                                    :: isub,j,j2,i,i2
         INTEGER,DIMENSION(5)                       :: ldl5,ldu5 
         INTEGER,DIMENSION(4)                       :: ldl4,ldu4 
         INTEGER,DIMENSION(4)                       :: ldl3,ldu3 
         INTEGER                                    :: iopt,topoid
         INTEGER                                    :: aa,bb,cc,dd,ee,gg,iface
#if __MESH_DIM == __3D
         INTEGER                                    :: k,k2
#endif
         INTEGER                                    :: mlevp1,ilda
         REAL(MK)                                   :: t0

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
         REAL(MK),DIMENSION(:,:),POINTER :: tuc
         REAL(MK),DIMENSION(:,:),POINTER :: puc
#elif __MESH_DIM == __3D
        REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
        REAL(MK),DIMENSION(:,:,:),POINTER :: puc
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
       REAL(MK),DIMENSION(:,:,:),POINTER :: tuc
       REAL(MK),DIMENSION(:,:,:),POINTER :: puc
#elif __MESH_DIM == __3D
       REAL(MK),DIMENSION(:,:,:,:),POINTER :: tuc
       REAL(MK),DIMENSION(:,:,:,:),POINTER :: puc
#endif
#endif

         !----------------------------------------------------------------------
         !Externals
         !----------------------------------------------------------------------

         !----------------------------------------------------------------------
         !Initialize
         !----------------------------------------------------------------------

         CALL substart('ppm_mg_prolong',t0,info)


         !----------------------------------------------------------------------
         !  Check arguments
         !----------------------------------------------------------------------
         IF (ppm_debug .GT. 0) THEN
             IF (mlev.LT.1) THEN
                   info = ppm_error_error
                   CALL ppm_error(ppm_err_argument,'ppm_mg_prolong',  &
      &                'level must be >0',__LINE__,info)
                   GOTO 9999
             ENDIF
         ENDIF
         !----------------------------------------------------------------------
         !Definition of necessary variables and allocation of arrays
         !----------------------------------------------------------------------
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




         !----------------------------------------------------------------------
         !Implementation
         !----------------------------------------------------------------------
         mlevp1 = mlev + 1


         IF (l_print) THEN
          WRITE(cbuf,*) 'WELCOME TO THE PROLONG LEVEL:',mlev
          CALL PPM_WRITE(ppm_rank,'mg_prolong',cbuf,info)
         ENDIF


#if __DIM == __SFIELD
#if __MESH_DIM == __2D           
         !----------------------------------------------------------------------
         !prolongation using a 9-point operator
         !----------------------------------------------------------------------


          DO isub=1,nsubs

          tuc=>mgfield(isub,mlevp1)%uc   

          puc=>mgfield(isub,mlev)%uc   


            DO j=1,max_node(2,mlevp1)
               j2=2*j
               DO i=1,max_node(1,mlevp1)
                  i2=2*i
                     puc(i2-1,j2-1) = &
      &                              puc(i2-1,j2-1) + &
      &                              tuc(i,j)

                     puc(i2,j2-1) = & 
      &                              puc(i2,j2-1) + &
      &                    0.5_MK * (tuc(i,j)+&
      &                              tuc(i+1,j))

                     puc(i2-1,j2) = &
      &                              puc(i2-1,j2) + &
      &                   0.5_MK * ( tuc(i,j) + & 
      &                              tuc(i,j+1))

                     puc(i2,j2)  = &
      &                             puc(i2,j2) + &
      &                  0.25_MK * ( tuc(i,j)+&
      &                              tuc(i+1,j) + &
      &                              tuc(i+1,j+1)+&
      &                              tuc(i,j+1))

               ENDDO
            ENDDO
         ENDDO


#elif __MESH_DIM == __3D

         !----------------------------------------------------------------------
         !prolongation using a 27-point operator
         !----------------------------------------------------------------------


          DO isub=1,nsubs
            tuc=>mgfield(isub,mlevp1)%uc 
            puc=>mgfield(isub,mlev)%uc 
            DO k=1,max_node(3,mlevp1)
               k2=2*k
               DO j=1,max_node(2,mlevp1)
                  j2=2*j
                  DO i=1,max_node(1,mlevp1)
                     i2=2*i
                        puc(i2-1,j2-1,k2-1) = &
      &                          puc(i2-1,j2-1,k2-1) + &
      &                          tuc(i,j,k)

                        puc(i2,j2-1,k2-1) = & 
      &                             puc(i2,j2-1,k2-1) + &
      &                   0.5_MK * (tuc(i,j,k)+&
      &                             tuc(i+1,j,k))

                        puc(i2-1,j2,k2-1) = &
      &                             puc(i2-1,j2,k2-1) + &
      &                  0.5_MK * ( tuc(i,j,k) + & 
      &                             tuc(i,j+1,k))

                        puc(i2,j2,k2-1)  = &
      &                            puc(i2,j2,k2-1) + &
      &                 0.25_MK * ( tuc(i,j,k)+&
      &                             tuc(i+1,j,k) + &
      &                             tuc(i+1,j+1,k)+&
      &                             tuc(i,j+1,k))


                        puc(i2-1,j2-1,k2) = &
      &                        puc(i2-1,j2-1,k2) + &
      &                 0.5_MK * (tuc(i,j,k)+&
      &                          tuc(i,j,k+1))



                        puc(i2,j2-1,k2)  = &
      &                             puc(i2,j2-1,k2) + &
      &                 0.25_MK * ( tuc(i,j,k)+&
      &                             tuc(i+1,j,k) + &
      &                             tuc(i+1,j,k+1)+&
      &                             tuc(i,j,k+1))

                        puc(i2-1,j2,k2)  = &
      &                             puc(i2-1,j2,k2) + &
      &                 0.25_MK * ( tuc(i,j,k)+&
      &                             tuc(i,j+1,k) + &
      &                             tuc(i,j,k+1)+&
      &                             tuc(i,j+1,k+1))


                       puc(i2,j2,k2)  = &
      &                             puc(i2,j2,k2) + &
      &                 0.125_MK * (tuc(i,j,k)+&
      &                            tuc(i+1,j,k) + &
      &                            tuc(i+1,j+1,k)+&
      &                            tuc(i,j+1,k)+&
      &                            tuc(i,j,k+1) + &
      &                            tuc(i+1,j,k+1)+&
      &                            tuc(i,j+1,k+1)+&
      &                            tuc(i+1,j+1,k+1))



                  ENDDO
               ENDDO
            ENDDO
         ENDDO


#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D           
         !----------------------------------------------------------------------
         !prolongation using a 9-point operator
         !----------------------------------------------------------------------


          DO isub=1,nsubs
             tuc=>mgfield(isub,mlevp1)%uc 
             puc=>mgfield(isub,mlev)%uc 

            DO j=1,max_node(2,mlevp1)
               j2=2*j
               DO i=1,max_node(1,mlevp1)
                  i2=2*i
                DO ilda=1,vecdim
                     puc(ilda,i2-1,j2-1) = &
      &                              puc(ilda,i2-1,j2-1) + &
      &                              tuc(ilda,i,j)

                     puc(ilda,i2,j2-1) = & 
      &                              puc(ilda,i2,j2-1) + &
      &                    0.5_MK * (tuc(ilda,i,j)+&
      &                              tuc(ilda,i+1,j))

                     puc(ilda,i2-1,j2) = &
      &                              puc(ilda,i2-1,j2) + &
      &                   0.5_MK * ( tuc(ilda,i,j) + & 
      &                              tuc(ilda,i,j+1))

                     puc(ilda,i2,j2)  = &
      &                             puc(ilda,i2,j2) + &
      &                  0.25_MK * ( tuc(ilda,i,j)+&
      &                              tuc(ilda,i+1,j) + &
      &                              tuc(ilda,i+1,j+1)+&
      &                              tuc(ilda,i,j+1))

                ENDDO 
               ENDDO
            ENDDO
         ENDDO


#elif __MESH_DIM == __3D

         !----------------------------------------------------------------------
         !prolongation using a 27-point operator
         !----------------------------------------------------------------------

          DO isub=1,nsubs
             tuc=>mgfield(isub,mlevp1)%uc 
             puc=>mgfield(isub,mlev)%uc 

            DO k=1,max_node(3,mlevp1)
               k2=2*k
               DO j=1,max_node(2,mlevp1)
                  j2=2*j
                  DO i=1,max_node(1,mlevp1)
                     i2=2*i
#ifdef __VECTOR

                        puc(1,i2-1,j2-1,k2-1) = &
      &                          puc(1,i2-1,j2-1,k2-1) + &
      &                          tuc(1,i,j,k)

                        puc(1,i2,j2-1,k2-1) = & 
      &                             puc(1,i2,j2-1,k2-1) + &
      &                   0.5_MK * (tuc(1,i,j,k)+&
      &                             tuc(1,i+1,j,k))

                        puc(1,i2-1,j2,k2-1) = &
      &                             puc(1,i2-1,j2,k2-1) + &
      &                  0.5_MK * ( tuc(1,i,j,k) + & 
      &                             tuc(1,i,j+1,k))

                        puc(1,i2,j2,k2-1)  = &
      &                            puc(1,i2,j2,k2-1) + &
      &                 0.25_MK * ( tuc(1,i,j,k)+&
      &                             tuc(1,i+1,j,k) + &
      &                             tuc(1,i+1,j+1,k)+&
      &                             tuc(1,i,j+1,k))


                        puc(1,i2-1,j2-1,k2) = &
      &                        puc(1,i2-1,j2-1,k2) + &
      &                 0.5_MK * (tuc(1,i,j,k)+&
      &                          tuc(1,i,j,k+1))



                        puc(1,i2,j2-1,k2)  = &
      &                             puc(1,i2,j2-1,k2) + &
      &                 0.25_MK * ( tuc(1,i,j,k)+&
      &                             tuc(1,i+1,j,k) + &
      &                             tuc(1,i+1,j,k+1)+&
      &                             tuc(1,i,j,k+1))

                        puc(1,i2-1,j2,k2)  = &
      &                             puc(1,i2-1,j2,k2) + &
      &                 0.25_MK * ( tuc(1,i,j,k)+&
      &                             tuc(1,i,j+1,k) + &
      &                             tuc(1,i,j,k+1)+&
      &                             tuc(1,i,j+1,k+1))


                       puc(1,i2,j2,k2)  = &
      &                             puc(1,i2,j2,k2) + &
      &                 0.125_MK * (tuc(1,i,j,k)+&
      &                            tuc(1,i+1,j,k) + &
      &                            tuc(1,i+1,j+1,k)+&
      &                            tuc(1,i,j+1,k)+&
      &                            tuc(1,i,j,k+1) + &
      &                            tuc(1,i+1,j,k+1)+&
      &                            tuc(1,i,j+1,k+1)+&
      &                            tuc(1,i+1,j+1,k+1))

                        puc(2,i2-1,j2-1,k2-1) = &
      &                          puc(2,i2-1,j2-1,k2-1) + &
      &                          tuc(2,i,j,k)

                        puc(2,i2,j2-1,k2-1) = & 
      &                             puc(2,i2,j2-1,k2-1) + &
      &                   0.5_MK * (tuc(2,i,j,k)+&
      &                             tuc(2,i+1,j,k))

                        puc(2,i2-1,j2,k2-1) = &
      &                             puc(2,i2-1,j2,k2-1) + &
      &                  0.5_MK * ( tuc(2,i,j,k) + & 
      &                             tuc(2,i,j+1,k))

                        puc(2,i2,j2,k2-1)  = &
      &                            puc(2,i2,j2,k2-1) + &
      &                 0.25_MK * ( tuc(2,i,j,k)+&
      &                             tuc(2,i+1,j,k) + &
      &                             tuc(2,i+1,j+1,k)+&
      &                             tuc(2,i,j+1,k))


                        puc(2,i2-1,j2-1,k2) = &
      &                        puc(2,i2-1,j2-1,k2) + &
      &                 0.5_MK * (tuc(2,i,j,k)+&
      &                          tuc(2,i,j,k+1))



                        puc(2,i2,j2-1,k2)  = &
      &                             puc(2,i2,j2-1,k2) + &
      &                 0.25_MK * ( tuc(2,i,j,k)+&
      &                             tuc(2,i+1,j,k) + &
      &                             tuc(2,i+1,j,k+1)+&
      &                             tuc(2,i,j,k+1))

                        puc(2,i2-1,j2,k2)  = &
      &                             puc(2,i2-1,j2,k2) + &
      &                 0.25_MK * ( tuc(2,i,j,k)+&
      &                             tuc(2,i,j+1,k) + &
      &                             tuc(2,i,j,k+1)+&
      &                             tuc(2,i,j+1,k+1))


                       puc(2,i2,j2,k2)  = &
      &                             puc(2,i2,j2,k2) + &
      &                 0.125_MK * (tuc(2,i,j,k)+&
      &                            tuc(2,i+1,j,k) + &
      &                            tuc(2,i+1,j+1,k)+&
      &                            tuc(2,i,j+1,k)+&
      &                            tuc(2,i,j,k+1) + &
      &                            tuc(2,i+1,j,k+1)+&
      &                            tuc(2,i,j+1,k+1)+&
      &                            tuc(2,i+1,j+1,k+1))


                        puc(3,i2-1,j2-1,k2-1) = &
      &                          puc(3,i2-1,j2-1,k2-1) + &
      &                          tuc(3,i,j,k)

                        puc(3,i2,j2-1,k2-1) = & 
      &                             puc(3,i2,j2-1,k2-1) + &
      &                   0.5_MK * (tuc(3,i,j,k)+&
      &                             tuc(3,i+1,j,k))

                        puc(3,i2-1,j2,k2-1) = &
      &                             puc(3,i2-1,j2,k2-1) + &
      &                  0.5_MK * ( tuc(3,i,j,k) + & 
      &                             tuc(3,i,j+1,k))

                        puc(3,i2,j2,k2-1)  = &
      &                            puc(3,i2,j2,k2-1) + &
      &                 0.25_MK * ( tuc(3,i,j,k)+&
      &                             tuc(3,i+1,j,k) + &
      &                             tuc(3,i+1,j+1,k)+&
      &                             tuc(3,i,j+1,k))


                        puc(3,i2-1,j2-1,k2) = &
      &                        puc(3,i2-1,j2-1,k2) + &
      &                 0.5_MK * (tuc(3,i,j,k)+&
      &                          tuc(3,i,j,k+1))



                        puc(3,i2,j2-1,k2)  = &
      &                             puc(3,i2,j2-1,k2) + &
      &                 0.25_MK * ( tuc(3,i,j,k)+&
      &                             tuc(3,i+1,j,k) + &
      &                             tuc(3,i+1,j,k+1)+&
      &                             tuc(3,i,j,k+1))

                        puc(3,i2-1,j2,k2)  = &
      &                             puc(3,i2-1,j2,k2) + &
      &                 0.25_MK * ( tuc(3,i,j,k)+&
      &                             tuc(3,i,j+1,k) + &
      &                             tuc(3,i,j,k+1)+&
      &                             tuc(3,i,j+1,k+1))


                       puc(3,i2,j2,k2)  = &
      &                             puc(3,i2,j2,k2) + &
      &                 0.125_MK * (tuc(3,i,j,k)+&
      &                            tuc(3,i+1,j,k) + &
      &                            tuc(3,i+1,j+1,k)+&
      &                            tuc(3,i,j+1,k)+&
      &                            tuc(3,i,j,k+1) + &
      &                            tuc(3,i+1,j,k+1)+&
      &                            tuc(3,i,j+1,k+1)+&
      &                            tuc(3,i+1,j+1,k+1))
#else
                   DO ilda=1,vecdim
                        puc(ilda,i2-1,j2-1,k2-1) = &
      &                          puc(ilda,i2-1,j2-1,k2-1) + &
      &                          tuc(ilda,i,j,k)

                        puc(ilda,i2,j2-1,k2-1) = & 
      &                             puc(ilda,i2,j2-1,k2-1) + &
      &                   0.5_MK * (tuc(ilda,i,j,k)+&
      &                             tuc(ilda,i+1,j,k))

                        puc(ilda,i2-1,j2,k2-1) = &
      &                             puc(ilda,i2-1,j2,k2-1) + &
      &                  0.5_MK * ( tuc(ilda,i,j,k) + & 
      &                             tuc(ilda,i,j+1,k))

                        puc(ilda,i2,j2,k2-1)  = &
      &                            puc(ilda,i2,j2,k2-1) + &
      &                 0.25_MK * ( tuc(ilda,i,j,k)+&
      &                             tuc(ilda,i+1,j,k) + &
      &                             tuc(ilda,i+1,j+1,k)+&
      &                             tuc(ilda,i,j+1,k))


                        puc(ilda,i2-1,j2-1,k2) = &
      &                        puc(ilda,i2-1,j2-1,k2) + &
      &                 0.5_MK * (tuc(ilda,i,j,k)+&
      &                          tuc(ilda,i,j,k+1))



                        puc(ilda,i2,j2-1,k2)  = &
      &                             puc(ilda,i2,j2-1,k2) + &
      &                 0.25_MK * ( tuc(ilda,i,j,k)+&
      &                             tuc(ilda,i+1,j,k) + &
      &                             tuc(ilda,i+1,j,k+1)+&
      &                             tuc(ilda,i,j,k+1))

                        puc(ilda,i2-1,j2,k2)  = &
      &                             puc(ilda,i2-1,j2,k2) + &
      &                 0.25_MK * ( tuc(ilda,i,j,k)+&
      &                             tuc(ilda,i,j+1,k) + &
      &                             tuc(ilda,i,j,k+1)+&
      &                             tuc(ilda,i,j+1,k+1))


                       puc(ilda,i2,j2,k2)  = &
      &                             puc(ilda,i2,j2,k2) + &
      &                 0.125_MK * (tuc(ilda,i,j,k)+&
      &                            tuc(ilda,i+1,j,k) + &
      &                            tuc(ilda,i+1,j+1,k)+&
      &                            tuc(ilda,i,j+1,k)+&
      &                            tuc(ilda,i,j,k+1) + &
      &                            tuc(ilda,i+1,j,k+1)+&
      &                            tuc(ilda,i,j+1,k+1)+&
      &                            tuc(ilda,i+1,j+1,k+1))

                    ENDDO
#endif
                  ENDDO
               ENDDO
            ENDDO
         ENDDO


#endif
#endif

         !----------------------------------------------------------------------
         !RETURN
         !----------------------------------------------------------------------

 9999    CONTINUE
         CALL substop('ppm_mg_prolong',t0,info)
         RETURN
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
       END SUBROUTINE ppm_mg_prolong_3d_vec_d
#endif
#endif
#endif


