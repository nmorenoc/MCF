#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_util_gmres
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Solve Av=b
      !
      !  Input        : A          (F) function that computes the LHS A.v
      !                 b          (F) Right-hand side
      !                 n          (I) size
      !                 m          (I) restart size
      !  
      !  Input/Output : X          (F) Initial guess/solution
      !                 tol        (F) tolerance
      !
      !  Output       : info       (I) return status.
      !
      !  Remarks      : 1) There is an optimal m...  No real rule to get it...
      !                    I would advise experimentation with values of m around
      !                    5 -> 10 -> 30, depending on the size of the problem
      !                    and the available memory
      !
      !  References   : Y. SAAD AND M. SCHULTZ, GMRES: A generalized minimal 
      !                 residual algorithm for solving nonsymmetric linear 
      !                 systems, SIAM J. Sci. Statist. Comput., 7 (1986), 
      !                 pp. 856-869.
      !
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_util_gmres.f,v $
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
      SUBROUTINE ppm_util_gmres_s(A,b,X,n,m,tol,info,maxiter)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_util_gmres_d(A,b,X,n,m,tol,info,maxiter)
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
      USE ppm_module_write
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

#ifdef HAVE_MPI
      INCLUDE 'mpif.h'
#else
#include "fakempi.h"
#endif

      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
#if     __KIND == __SINGLE_PRECISION
      INTERFACE
         FUNCTION A(in,out,ninout,info)
         INTEGER                          , INTENT(OUT) :: info
         REAL(KIND(1.0E0)), DIMENSION(:),   INTENT(IN)  :: in
         REAL(KIND(1.0E0)), DIMENSION(:),   INTENT(OUT) :: out
         INTEGER                          , INTENT(IN)  :: ninout
         INTEGER                                        :: A
         END FUNCTION A
      END INTERFACE
#else 
      INTERFACE
         FUNCTION A(in,out,ninout,info)
         INTEGER                          , INTENT(OUT) :: info
         REAL(KIND(1.0D0)), DIMENSION(:),   INTENT(IN)  :: in
         REAL(KIND(1.0D0)), DIMENSION(:),   INTENT(OUT) :: out
         INTEGER                          , INTENT(IN)  :: ninout
         INTEGER                                        :: A
         END FUNCTION A
      END INTERFACE
#endif
      REAL(MK), DIMENSION(:), POINTER       :: x
      REAL(MK), DIMENSION(:), POINTER       :: b
      REAL(MK),               INTENT(INOUT) :: tol
      INTEGER               , INTENT(IN)    :: n
      INTEGER               , INTENT(IN)    :: m
      INTEGER               , INTENT(OUT)   :: info
      INTEGER, OPTIONAL     , INTENT(IN)    :: maxiter
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(MK)                              :: t0,tolin,lmyeps
      INTEGER                               :: solver_info
      
      ! Accumulators
      ! Need double precision for the inner products
      REAL(KIND(1.0D0)),DIMENSION(:)  ,POINTER :: tmphd
      REAL(KIND(1.0D0)),DIMENSION(:,:),POINTER :: ghd
      REAL(KIND(1.0D0))                        :: lnorm,gnorm,groti,grotip1
      REAL(KIND(1.0D0)),DIMENSION(:),POINTER   :: s,y
      REAL(KIND(1.0D0)),DIMENSION(:,:),POINTER :: cossin
      
      REAL(MK), DIMENSION(:), POINTER       :: w
      REAL(MK), DIMENSION(:,:), POINTER     :: v
      INTEGER                               :: i,j,k,l,throwaway, iopt
      INTEGER, DIMENSION(3)                 :: ldl,ldu
      LOGICAL                               :: breakdown
      CHARACTER(LEN=ppm_char)               :: mesg
#ifdef USE_MPI
      INTEGER                               :: MPTYPE
#endif


#ifdef USE_MPI
      !-------------------------------------------------------------------------
      !  Define MPI data type
      !-------------------------------------------------------------------------
      MPTYPE = MPI_DOUBLE_PRECISION
#endif

      !-----------------------------------------------------------------------
      !  call substart
      !-----------------------------------------------------------------------
      CALL substart('ppm_util_gmres',t0,info)

      !-----------------------------------------------------------------------
      !  check input arguments
      !-----------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         IF (n.LT.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_gmres',  &
     &            'passed a negative size n',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      IF (ppm_debug .GT. 0) THEN
         IF (m.LT.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_gmres',  &
     &            'passed a restart parameter m < 1',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      IF (ppm_debug .GT. 0) THEN
         IF (m.GT.n) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_gmres',  &
     &            'restart parameter m > size n',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Allocation
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldl(1) = 1
      ldu(1) = n
      CALL ppm_alloc(w,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'allocation of vector w',__LINE__,info)
          GOTO 9999
      ENDIF
      
      iopt = ppm_param_alloc_fit
      ldl(1) = 1
      ldu(1) = m+1
      CALL ppm_alloc(s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'allocation of vector s',__LINE__,info)
          GOTO 9999
      ENDIF
      
      iopt = ppm_param_alloc_fit
      ldl(1) = 1
      ldu(1) = n
      ldl(2) = 1
      ldu(2) = m+1
      CALL ppm_alloc(v,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'allocation of basis vectors v',__LINE__,info)
          GOTO 9999
      ENDIF
      iopt = ppm_param_alloc_fit
      ldl(1) = 1
      ldu(1) = m
      CALL ppm_alloc(y,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'allocation of vectors y',__LINE__,info)
          GOTO 9999
      ENDIF

      iopt = ppm_param_alloc_fit
      ldl(1) = 1
      ldu(1) = m+1
      CALL ppm_alloc(tmphd,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'allocation of vector TMPH of local inner products',__LINE__,info)
          GOTO 9999
      ENDIF
      
      iopt = ppm_param_alloc_fit
      ldl(1) = 1
      ldu(1) = m+1
      ldl(2) = 1
      ldu(2) = m
      CALL ppm_alloc(ghd,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'allocation of matrix GH of global inner products',__LINE__,info)
          GOTO 9999
      ENDIF

      iopt = ppm_param_alloc_fit
      ldl(1) = 1
      ldu(1) = 2
      ldl(2) = 1
      ldu(2) = m
      CALL ppm_alloc(cossin,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'allocation of matrix COSSIN of cosines/sines of rotation angles',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Iterations
      !-------------------------------------------------------------------------
      tolin = tol
      j = 0
      DO
         j = j+1
         !----------------------------------------------------------------------
         !  (Re)-Initialization
         !----------------------------------------------------------------------
         throwaway = A(x,w,n,info)
         w = b(1:n) - w(1:n)
#ifdef USE_MPI
         lnorm = dot_product(w,w)
         !lnorm = 0.0
         !DO k=1,n
         !   lnorm = lnorm +  REAL(w(k),KIND(1.0D0))**2
         !END DO
         CALL MPI_Allreduce(lnorm,gnorm,1,MPTYPE,MPI_SUM,ppm_comm,info)
#else
         gnorm = dot_product(w,w)
         !gnorm = 0.0
         !DO k=1,n
         !   gnorm = gnorm + REAL(w(k),KIND(1.0D0))**2
         !END DO
#endif
         gnorm = SQRT(gnorm)
         v(:,1) = w / gnorm
         
         s = 0.0
         s(1) = gnorm
         
         breakdown = .FALSE.
         ghd = 0.0
         DO i=1,m
            !-------------------------------------------------------------------
            !  Growing the Krylov space basis: Modified Gram-Schmidt
            !  NB: could switch to a more stable orthognalization scheme
            !      like Householder 
            !-------------------------------------------------------------------
            throwaway = A(v(:,i),w,n,info)
#ifdef USE_MPI
            tmphd = 0.0
            DO  k=1,i
               tmphd(k) = dot_product(w,v(:,k))
               !DO l=1,n
               !   tmphd(k) = tmphd(k) + &
    !&                        REAL(w(l),KIND(1.0D0))*REAL(v(l,k),KIND(1.0D0))
               !END DO
            END DO
            CALL MPI_Allreduce(tmphd(1:i),ghd(1:i,i),i,MPTYPE,MPI_SUM,ppm_comm,info)
#else
            DO  k=1,i
               ghd(k,i) = dot_product(w,v(:,k))
               !DO l=1,n
               !   ghd(k,i) = ghd(k,i) + &
    !&                        REAL(w(l),KIND(1.0D0))*REAL(v(l,k),KIND(1.0D0))
               !END DO
            END DO
#endif
            DO  k=1,i
               w = w - ghd(k,i)*v(:,k)
            END DO

#ifdef USE_MPI
            lnorm = dot_product(w,w)
            !lnorm = 0.0
            !DO k=1,n
            !   lnorm = lnorm + REAL(w(k),KIND(1.0D0))**2
            !END DO
            CALL MPI_Allreduce(lnorm,gnorm,1,MPTYPE,MPI_SUM,ppm_comm,info)
#else
            gnorm = dot_product(w,w)
            !gnorm = 0.0
            !DO k=1,n
            !   gnorm = gnorm + REAL(w(k),KIND(1.0D0))**2
            !END DO
#endif
            gnorm = SQRT(gnorm)
            ghd(i+1,i) = gnorm

            IF (ghd(i+1,i).LT.lmyeps) THEN
               IF (ppm_debug .GE. 1) THEN
                  WRITE(mesg,'(A, ES9.2,A,ES9.2)') 'Breakdown: ', ghd(i+1,i), ' < ', lmyeps
                  CALL ppm_write(ppm_rank,'ppm_util_gmres',  &
                  &              mesg,info)
               ENDIF
               breakdown = .TRUE.
               EXIT
            END IF

            v(:,i+1) = w / ghd(i+1,i)
            
            ! Apply previous Givens rotations 
            DO k=1,i-1
               groti   =  cossin(1,k)*ghd(k,i) + cossin(2,k)*ghd(k+1,i)
               grotip1 = -cossin(2,k)*ghd(k,i) + cossin(1,k)*ghd(k+1,i)
               ghd(k,i) = groti
               ghd(k+1,i) = grotip1
            END DO
            !WRITE(*,*) 'gh(',i,'+1,',i,') after prev rots= ', gh(i+1,i)
            ! Build the new one
            lnorm = SQRT(ghd(i,i)**2+ghd(i+1,i)**2)
            cossin(1,i)=ghd(i,i)  /lnorm
            cossin(2,i)=ghd(i+1,i)/lnorm
            ! And apply it
            groti   =  cossin(1,i)*ghd(i,i) + cossin(2,i)*ghd(i+1,i)
            grotip1 = -cossin(2,i)*ghd(i,i) + cossin(1,i)*ghd(i+1,i)
            ghd(i,i) = groti
            ghd(i+1,i) = 0.0
            ! Update s
            groti   =  cossin(1,i)*s(i)! + cossin(2,i)*s(i+1)
            grotip1 = -cossin(2,i)*s(i)! + cossin(1,i)*s(i+1)
            s(i) = groti
            s(i+1) = grotip1
            tol = ABS(s(i+1))
            
            IF (tol.LT.tolin) THEN
               ! We are done, update tol, compute solution and exit
               IF (ppm_debug .GE. 1) THEN
                  WRITE(mesg,'(A, ES9.2,A,ES9.2,A,I4)') 'Converged: ', &
    &                   tol, ' < ', tolin, ', Outer iterations: ',j
                  CALL ppm_write(ppm_rank,'ppm_util_gmres',  &
                  &              mesg,info)
               ENDIF
               
               CALL ppm_util_gmres_solveupper(ghd,s,y,i,info)
               IF (info .NE. 0) THEN
                  info = ppm_error_fatal
                  CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',       &
     &                 'failed to solve hy=s',__LINE__,info)
                  GOTO 9999
               ENDIF

               DO k=1,i
                  x = x + y(k)*v(1:n,k)
               END DO
               solver_info = ppm_gmres_param_success
               GOTO 8000
            END IF
            
         END DO
         
         ! We got out of the inner loop, either i=m+1 or breakdown and i<=m+1
         tol = ABS(s(i))
         IF (.NOT.breakdown) THEN
            i = m
         END IF
         CALL ppm_util_gmres_solveupper(ghd,s,y,i,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_div_zero,'ppm_util_gmres',        &
     &           'failed to solve hy=s',__LINE__,info)
            solver_info = ppm_gmres_param_failure
            GOTO 8000
         ENDIF
         DO k=1,i
            x = x + y(k)*v(1:n,k)
         END DO
         IF (tol.LT.tolin) THEN
            IF (ppm_debug .GE. 1) THEN
               WRITE(mesg,'(A, ES9.2,A,ES9.2,A,I4)') 'Converged: ', &
     &               tol, ' < ', tolin,', Outer iterations: ',j
               CALL ppm_write(ppm_rank,'ppm_util_gmres',  &
               &              mesg,info)
            ENDIF
            solver_info = ppm_gmres_param_success
            GOTO 8000
         END IF
         IF (PRESENT(maxiter)) THEN
            IF (j.GE.maxiter) THEN
               IF (ppm_debug .GE. 1) THEN
                  WRITE(mesg,'(A,I4,A,ES9.2,A,ES9.2)') 'Reached ', &
     &                 maxiter,' iterations and residue still at ', &
     &                 tol, ' < ', tolin
                  CALL ppm_write(ppm_rank,'ppm_util_gmres',  &
                  &              mesg,info)
               ENDIF
               solver_info = ppm_gmres_param_maxiter
               GOTO 8000
            END IF
         ENDIF
      END DO

      
      !-------------------------------------------------------------------------
      !  De-allocate
      !-------------------------------------------------------------------------
 8000 CONTINUE
      iopt = ppm_param_dealloc
      CALL ppm_alloc(cossin,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'deallocation of vector tmpv',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'deallocation of vector s',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(v,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'deallocation of basis vectors v',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(w,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'deallocation of basis vectors v',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(y,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'deallocation of vectors y',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(tmphd,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'deallocation of vector TMPH of local inner products',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ghd,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_util_gmres',        &
     &        'deallocation of matrix GH of global inner products',__LINE__,info)
          GOTO 9999
      ENDIF
      !-------------------------------------------------------------------------
      !  Made it through: we set the output info to the status of the solver
      !-------------------------------------------------------------------------
      info = solver_info

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_gmres',t0,info)
      RETURN

#if   __KIND == __SINGLE_PRECISION
END SUBROUTINE ppm_util_gmres_s
#elif __KIND == __DOUBLE_PRECISION
END SUBROUTINE ppm_util_gmres_d
#endif
