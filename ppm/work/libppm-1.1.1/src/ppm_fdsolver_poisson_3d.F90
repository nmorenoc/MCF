#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_poisson_3d.f
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine solves the poisson equation with the data 
      !               as right hand side in fourier space
      !                  - Laplacian of Phi = omega 
      !                 
      !
      !  Input        :  
      !                 lda(3)          (I)      size of local data field
      !                 istart(3)       (I)      starting index of local field
      !                 Nm(3)           (I)      size of global data
      !                 length(3)       (I)      length of the domain 
      !                                
      !
      !  Input/output : 
      !                 fdata(:,:)       (F)     data in fourier space
      !                                         right hand side  of poisson
      !                                         equation   
      !                                
      !
      !  Output       : 
      !                 info       (       I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_poisson_3d.f,v $
      !  Revision 1.6  2006/09/04 18:34:44  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.5  2004/10/01 16:08:58  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.4  2004/09/24 11:32:20  hiebers
      !  added MK to float numbers, exchanged data by fdata
      !
      !  Revision 1.3  2004/07/26 13:49:17  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.2  2004/07/26 11:59:38  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.1  2004/07/26 08:52:28  hiebers
      !  Recommited, formerly ppm_module_fieldsolver
      !
      !  Revision 1.1  2004/05/19 15:36:00  hiebers
      !  implementation from scratch
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if __KIND == __COMPLEX
      SUBROUTINE ppm_fdsolver_poisson_3dc(fdata,lda,istart,length,Nm,info)
#elif __KIND == __DOUBLE_COMPLEX
      SUBROUTINE ppm_fdsolver_poisson_3dcc(fdata,lda,istart,length,Nm,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_write
      USE ppm_module_error

      IMPLICIT NONE

#if   __KIND == __COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif



      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data
      COMPLEX(MK), DIMENSION(:,:,:),    POINTER          :: fdata
      INTEGER, DIMENSION(3),            INTENT(IN)       :: lda, istart,Nm
      REAL(MK), DIMENSION(3),           INTENT(IN)       :: length
      INTEGER,                          INTENT(  OUT)    :: info


      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                :: t0
      ! counters
      INTEGER                                 :: i,j,k
      INTEGER                                 :: i_global,j_global, k_global
      !Size of the data 
      INTEGER, DIMENSION(2)                   :: lda_end
      INTEGER                                 :: lda_end2
      ! wave number
      REAL(MK)                                :: kx, ky, kz
      REAL(MK)                                :: pi2_Lx, pi2_Ly, pi2_Lz
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_poisson_3d',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------

      IF (ppm_debug .GT. 0) THEN
          IF (lda(1).LE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_poisson_3d',  &
     &            ' mesh size: Nx must be >0 ' ,__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(2).LE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_poisson_3d',  &
     &            ' mesh size: Ny must be >0 ' ,__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(3).LE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_poisson_3d',  &
     &            ' mesh size: Ny must be >0 ' ,__LINE__,info)
              GOTO 9999
          ENDIF
       ENDIF


      !-------------------------------------------------------------------------
      !  Solve Poisson Equation
      !-------------------------------------------------------------------------

#if   __KIND == __COMPLEX 
      pi2_Lx = 2.0_MK*ppm_pi_s/length(1)
      pi2_Ly = 2.0_MK*ppm_pi_s/length(2)
      pi2_Lz = 2.0_MK*ppm_pi_s/length(3)
#elif __KIND == __DOUBLE_COMPLEX 
      pi2_Lx = 2.0_MK*ppm_pi_d/length(1)
      pi2_Ly = 2.0_MK*ppm_pi_d/length(2)
      pi2_Lz = 2.0_MK*ppm_pi_d/length(3)
#endif

      lda_end(1)=Nm(1)/2 +1
      lda_end(2)=Nm(2)/2 +1
      

      !-------------------------------------------------------------------------
      !  loop over all elements in fourier space and multiply 1/(kx^2+ky^2+kz^2)
      !-------------------------------------------------------------------------

      DO k=1,lda(3)
        DO j=1,lda(2)
          DO i=1,lda(1)
            ! decremented global index 
            i_global = istart(1)+i-2
            j_global = istart(2)+j-2
            k_global = istart(3)+k-2
            IF((i_global.EQ.0) .AND. (j_global.EQ.0).AND. (k_global.EQ.0) ) THEN
               fdata(i,j,k)= 0.0_MK
            ELSE
               !  consider complex conjugate elements (i> N/2+1)
               IF( i_global.GE.lda_end(1)) THEN
               i_global = Nm(1)-i_global
               ENDIF

               !  consider complex conjugate elements (j> N/2+1)
               IF( j_global.GE.lda_end(2)) THEN
               j_global = Nm(2)-j_global
               ENDIF

               kx=pi2_Lx*real(i_global,MK)  
               ky=pi2_Ly*real(j_global,MK)
               kz=pi2_Lz*real(k_global,MK)
               fdata(i,j,k)=-1.0_MK /(kx*kx + ky*ky + kz*kz)*fdata(i,j,k)
            ENDIF

           ENDDO
         ENDDO
      ENDDO



      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fdsolver_poisson_3d',t0,info)

      RETURN
#if __KIND == __COMPLEX
      END SUBROUTINE ppm_fdsolver_poisson_3dc
#elif __KIND == __DOUBLE_COMPLEX
      END SUBROUTINE ppm_fdsolver_poisson_3dcc
#endif
