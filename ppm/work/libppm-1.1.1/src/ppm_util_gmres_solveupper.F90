#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_util_gmres_solveupper
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Solve Ux=b
      !
      !  Input        : U          (F) upper triangular matrix
      !                 b          (F) Right-hand side
      !                 n          (I) size
      !
      !  Output       : info       (I) return status.
      !                 x          (F) Solution
      !
      !  Remarks      :
      !
      !  References   : 
      !
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_gmres_solveupper.f,v $
      !  Revision 1.1  2006/05/11 10:27:00  pchatela
      !  Initial insertion
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_util_gmres_solveupper_s(U,b,x,n,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_gmres_solveupper_d(U,b,x,n,info)
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN)  :: U
      REAL(MK), DIMENSION(:),   POINTER     :: x
      REAL(MK), DIMENSION(:),   INTENT(IN)  :: b
      INTEGER               , INTENT(IN)    :: n

      INTEGER               , INTENT(OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                              :: t0, lmyeps
      REAL(KIND(1.0D0))                     :: sum
      INTEGER                               :: i,j

      !-----------------------------------------------------------------------
      !  call substart
      !-----------------------------------------------------------------------
      CALL substart('ppm_util_gmres_solveupper',t0,info)

      !-----------------------------------------------------------------------
      !  check input arguments
      !-----------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         IF (n.LT.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_gmres_solveupper',  &
     &            'passed a negative size n',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif
      
      !-----------------------------------------------------------------------
      !  Elimination
      !-----------------------------------------------------------------------
      DO i=n,1,-1
         sum = b(i)
         DO j = i+1,n
            sum = sum - U(i,j)*x(j)
         END DO
         IF (ABS(U(i,i)).GT.lmyeps) THEN
            x(i) = sum/U(i,i)
         ELSE
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_util_gmres_solveupper',  &
     &            'a pivot of U is null',__LINE__,info)
            GOTO 9999
         END IF
      END DO
      
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_gmres_solveupper',t0,info)
      RETURN

#if   __KIND == __SINGLE_PRECISION
END SUBROUTINE ppm_util_gmres_solveupper_s
#elif __KIND == __DOUBLE_PRECISION
END SUBROUTINE ppm_util_gmres_solveupper_d
#endif
