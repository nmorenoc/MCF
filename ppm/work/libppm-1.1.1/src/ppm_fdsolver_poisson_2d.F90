#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_poisson_2d.f
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine solves the poisson equation with the data 
      !               as right hand side in fourier space
      !                  - Laplacian of Phi = omega 
      !                 
      !
      !  Input        :  
      !                 lda(2)          (I)      size of local data field
      !                 istart(2)       (I)      starting index of local field
      !                 Nm(2)           (I)      size of global data
      !                 length(2)       (I)      length of the domain 
      !                                
      !
      !  Input/output : 
      !                 fdat(:,:)       (F)     data in fourier space:
      !                                         right hand side  of poisson
      !                                         equation   
      !                                
      !
      !  Output       : 
      !                 info       (       I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_poisson_2d.f,v $
      !  Revision 1.7  2006/09/04 18:34:44  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.6  2005/02/16 12:02:50  hiebers
      !  simplified loop
      !
      !  Revision 1.5  2004/10/01 16:08:58  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.4  2004/07/26 15:38:46  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.3  2004/07/26 13:49:17  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.2  2004/07/26 11:59:38  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.1  2004/07/26 08:52:27  hiebers
      !  Recommited, formerly ppm_module_fieldsolver
      !
      !  Revision 1.1  2004/05/19 15:35:25  hiebers
      !  implementation from scratch
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if __KIND == __COMPLEX
      SUBROUTINE ppm_fdsolver_poisson_2dc(fdat, lda, istart, length, Nm, info)
#elif __KIND == __DOUBLE_COMPLEX
      SUBROUTINE ppm_fdsolver_poisson_2dcc(fdat, lda, istart, length, Nm, info)
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
      COMPLEX(MK), DIMENSION(:,:),      INTENT(INOUT)    :: fdat
      INTEGER, DIMENSION(2),            INTENT(IN)       :: lda, istart, Nm
      REAL(MK), DIMENSION(2),           INTENT(IN)       :: length
      INTEGER,                          INTENT(  OUT)    :: info


      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                :: t0
      ! counters
      INTEGER                                 :: i,j
      INTEGER                                 :: i_global,j_global
      !Size of the data 
      INTEGER                                 :: lda_end
      ! wave number
      REAL(MK)                                :: kx, ky
      REAL(MK)                                :: pi2_Lx, pi2_Ly
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_poisson_2d',t0,info)

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
       ENDIF


      !-------------------------------------------------------------------------
      !  Solve Poisson Equation
      !-------------------------------------------------------------------------

#if   __KIND == __COMPLEX 
      pi2_Lx = 2.0_MK*ppm_pi_s/length(1)
      pi2_Ly = 2.0_MK*ppm_pi_s/length(2)
#elif __KIND == __DOUBLE_COMPLEX 
      pi2_Lx = 2.0_MK*ppm_pi_d/length(1)
      pi2_Ly = 2.0_MK*ppm_pi_d/length(2)
#endif

      lda_end=Nm(1)/2 +1


      DO j=1,lda(2)
        DO i=1,lda(1)
            ! decremented global index 
            i_global = istart(1)+i-2
            j_global = istart(2)+j-2
            IF((i_global.EQ.0) .AND. (j_global.EQ.0) ) THEN
               fdat(i,j)= 0.0_MK
            ELSE
               !  consider complex conjugate elements (i> N/2+1)
               IF( i_global.GE.lda_end) THEN
               i_global = Nm(1)-i_global
               ENDIF

               kx=pi2_Lx*real(i_global,MK)  
               ky=pi2_Ly*real(j_global,MK)
               !print*, i_global, j_global, lda_end, Nm
               fdat(i,j)=-1.0_MK/(kx*kx + ky*ky)*fdat(i,j)
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fdsolver_poisson_2d',t0,info)

      RETURN
#if __KIND == __COMPLEX
      END SUBROUTINE ppm_fdsolver_poisson_2dc
#elif __KIND == __DOUBLE_COMPLEX
      END SUBROUTINE ppm_fdsolver_poisson_2dcc
#endif
