#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_step
      !-------------------------------------------------------------------------
      !
      !  Purpose      : integrates one stage of a mode
      !
      !  Input        : odeid                  (I) id of mode to advance
      !                 xp(:,:)                (F) xp (trad. coords)
      !                 up(:,:)                (F) function to integrate
      !                 dup(:,:)               (F) derivative of up (in time)
      !                 lda                    (I) leading dimension of up
      !                 Npart                  (I) numper of particles
      !                 bfr(:,:)               (F) buffer
      !                 istage                 (I) which stage were in
      !                 time(4)                (F) time data: time(1) is the
      !                                            start time, time(2) the end
      !                                            time, time(3) the current
      !                                            time and time (4) the
      !                                            time step
      !                 rhsfunc                (I) function pointer to the rhs
      !                 ipackdata(:,:)         (I) optional storage array
      !                 lpackdata(:,:)         (L) optional storage array
      !                 rpackdata(:,:)         (F) optional storage array
      !
      !  Output       : info                   (I) return status
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_ode_step.f,v $
      !  Revision 1.20  2006/09/04 18:34:54  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.19  2006/02/03 09:37:06  ivos
      !  Changed input arguments from INOUT to POINTER to allow for ghost_put.
      !
      !  Revision 1.18  2005/08/05 09:30:06  ivos
      !  Added the flabbergsting STS of Michael.
      !
      !  Revision 1.17  2005/07/22 08:03:43  pchatela
      !  Used a DO loop rather than a vector notation for the 1st step inside RK4.
      !  2 reasons:
      !  1) DO loops are used everywhere else in the file...
      !  2) It apparently causes a seg fault for a system with moderate memory (like my personal laptop)
      !
      !  Revision 1.16  2005/01/05 16:27:18  ivos
      !  fixed INTENT and POINTER attributes of the arguments.
      !
      !  Revision 1.15  2004/10/11 06:53:14  hiebers
      !  cosmetics
      !
      !  Revision 1.14  2004/10/01 15:15:11  hiebers
      !  added Runge Kutta 4th order
      !
      !  Revision 1.13  2004/08/26 15:21:00  michaebe
      !  inserted 1:npart to be able to handle ghosts.  ghosts?
      !
      !  Revision 1.12  2004/08/13 15:30:45  michaebe
      !  added mid point runge kutta
      !
      !  Revision 1.11  2004/08/12 13:45:23  michaebe
      !  modified ppm_ode_sent() = ..
      !
      !  Revision 1.10  2004/07/27 10:29:20  michaebe
      !  renamed the scheme parameters from ppm_ode_scheme.. to
      !  ppm_param_ode_scheme..
      !
      !  Revision 1.9  2004/07/26 14:57:18  michaebe
      !  renamed the vector and scalar defines
      !
      !  Revision 1.8  2004/07/26 13:49:19  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.7  2004/07/26 11:49:01  michaebe
      !  added overloaded end subroutine
      !
      !  Revision 1.6  2004/07/26 11:33:04  michaebe
      !  inserted the use of the ppm ode data module.
      !
      !  Revision 1.5  2004/07/26 08:13:51  michaebe
      !  Added scalar code for lda=1.
      !
      !  Revision 1.4  2004/07/26 07:59:51  michaebe
      !  Cosmetics. Lines partially too long.
      !
      !  Revision 1.3  2004/07/26 07:54:57  michaebe
      !  Atomized. Furthermore added carrier arrays to transport auxillary data
      !  to the user-given right hand side function.
      !
      !  Revision 1.2  2004/06/10 16:20:04  ivos
      !  Moved all xpp directtives to column 1. The NEC xpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.1  2004/02/19 08:33:57  michaebe
      !  initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __MODE == __SCA
      
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_ode_step_ss(odeid,xp,up,dup,lda,Npart,bfr,istage,time,&
      & rhsfunc, ipackdata, rpackdata, lpackdata, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_ode_step_ds(odeid,xp,up,dup,lda,Npart,bfr,istage,time,&
      & rhsfunc, ipackdata, rpackdata, lpackdata, info)
#endif

#elif __MODE == __VEC

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_ode_step_sv(odeid,xp,up,dup,lda,Npart,bfr,istage,time,&
      & rhsfunc, ipackdata, rpackdata, lpackdata, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_ode_step_dv(odeid,xp,up,dup,lda,Npart,bfr,istage,time,&
      & rhsfunc, ipackdata, rpackdata, lpackdata, info)
#endif

#endif
        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
#include "ppm_define.h"

        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_substart
        USE ppm_module_data_ode
        USE ppm_module_substop
        USE ppm_module_data
        USE ppm_module_error

        
        IMPLICIT NONE
#if     __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: mk = ppm_kind_single
#else
        INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
#if     __KIND == __SINGLE_PRECISION
        INTERFACE
           FUNCTION rhsfunc(xp,up,dup,lda,npart,ipack,&
                &lpack,rpack,info)
             INTEGER                          , INTENT(IN)  :: lda,npart
             INTEGER                          , INTENT(OUT) :: info
#if     __MODE == __SCA
             REAL(KIND(1.0E0)), DIMENSION(:,:), POINTER     :: xp
             REAL(KIND(1.0E0)), DIMENSION(:),   POINTER     :: up
             REAL(KIND(1.0E0)), DIMENSION(:),   POINTER     :: dup
#elif   __MODE == __VEC
             REAL(KIND(1.0E0)), DIMENSION(:,:), POINTER     :: xp,up
             REAL(KIND(1.0E0)), DIMENSION(:,:), POINTER     :: dup
#endif
      
             INTEGER,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: ipack
             LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lpack
             REAL(kind(1.0E0)), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rpack
             INTEGER                                     :: rhsfunc
           END FUNCTION rhsfunc
        END INTERFACE
#else
        INTERFACE
           FUNCTION rhsfunc(xp,up,dup,lda,npart,ipack,&
                &lpack,rpack,info)
             INTEGER                          , INTENT(IN)  :: lda,npart
             INTEGER                          , INTENT(OUT) :: info
#if     __MODE == __SCA
             REAL(KIND(1.0D0)), DIMENSION(:,:), POINTER     :: xp
             REAL(KIND(1.0D0)), DIMENSION(:),   POINTER     :: up
             REAL(KIND(1.0D0)), DIMENSION(:),   POINTER     :: dup
#elif   __MODE == __VEC
             REAL(KIND(1.0D0)), DIMENSION(:,:), POINTER     :: xp,up
             REAL(KIND(1.0D0)), DIMENSION(:,:), POINTER     :: dup
#endif
             INTEGER,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: ipack
             LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lpack
             REAL(kind(1.0D0)), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rpack
             INTEGER                                     :: rhsfunc
           END FUNCTION rhsfunc
        END INTERFACE
#endif
        

        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER,                    INTENT(  out) :: info
        INTEGER,                    INTENT(in   ) :: odeid
#if     __MODE == __SCA
        REAL(mk), DIMENSION(:  ), POINTER         :: up,dup        
#elif   __MODE == __VEC        
        REAL(mk), DIMENSION(:,:), POINTER         :: up,dup
#endif        
        REAL(mk), DIMENSION(:,:), POINTER         :: xp
        REAL(mk), DIMENSION(:,:), POINTER         :: bfr
        INTEGER,                    INTENT(in   ) :: istage
        REAL(mk), DIMENSION(4),     INTENT(inout) :: time

        INTEGER,                    INTENT(in   )   :: lda, Npart
        INTEGER,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: ipackdata
        REAL(mk), DIMENSION(:,:), INTENT(IN), OPTIONAL :: rpackdata
        LOGICAL,  DIMENSION(:,:), INTENT(IN), OPTIONAL :: lpackdata

        !-----------------------------------------------------------------------
        !  Local Variables
        !-----------------------------------------------------------------------
        REAL(mk)                                  :: t, dt
        INTEGER                                   :: scheme, throwaway, i, j
        INTEGER                                   :: umidmin, umidmax
        INTEGER                                   :: mid, ilda
        REAL(mk)                                  :: M_PI
        INTEGER                                   :: stsn
        REAL(mk), DIMENSION(20)                   :: stsnu
        REAL(mk)                                  :: tau
        !-----------------------------------------------------------------------
        !  fill the nu parameters for the sts scheme
        !-----------------------------------------------------------------------
        M_PI = ACOS(-1.0_MK)
        stsnu(1)  = 0.0_mk
        stsnu(5)  = 0.04_mk
        stsnu(7)  = 0.0015_mk
        stsnu(9)  = 0.04_mk
        stsnu(20) = 0.006_mk

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart('ppm_ode_step',t0,info)
        
        !-----------------------------------------------------------------------
        !  check input arguments
        !-----------------------------------------------------------------------
        IF(Npart.EQ.0) THEN
           !--------------------------------------------------------------------
           ! best case: nothing to do
           !--------------------------------------------------------------------
           time(3) = time(3) + time(4)
           GOTO 9999
        END IF
        
        IF(ppm_debug.GT.0) THEN
           IF(.NOT.ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_ode_step',&
                   & 'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
           END IF

           !--------------------------------------------------------------------
           ! check odeid
           !--------------------------------------------------------------------
           umidmin = LBOUND(ppm_internal_mid,1)
           umidmax = UBOUND(ppm_internal_mid,1)
           IF(odeid.LT.umidmin.OR.odeid.GT.umidmax) THEN
              !-----------------------------------------------------------------
              ! user mid does not exist
              !-----------------------------------------------------------------
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'ODEID does not exist',__LINE__,info)
              GOTO 9999
           ELSE
              IF(ppm_internal_mid(odeid).EQ.-HUGE(odeid)) THEN
                 !--------------------------------------------------------------
                 ! user mid does not exist
                 !--------------------------------------------------------------
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_ode_step',& 
                      & 'ODEID does not exist',__LINE__,info)
                 GOTO 9999
              END IF
           END IF
           
           
           !--------------------------------------------------------------------
           ! check dimension
           !--------------------------------------------------------------------
           IF(Npart.LT.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'Npart cannot be <0',__LINE__,info)
              GOTO 9999
           END IF
           IF(lda.LE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'LDA must be >00',__LINE__,info)
              GOTO 9999
           END IF
           IF(time(4).LE.0.0_mk) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'dt must be >=0',__LINE__,info)
              GOTO 9999
           END IF
           IF(time(3).LT.time(1)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'time must be >= tstart',__LINE__,info)
              GOTO 9999
           END IF

           IF(scheme.EQ.ppm_param_ode_scheme_sts) THEN
                IF(.NOT.PRESENT(ipackdata)) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                         & 'for STS you need to specify N in ipackdata(1,1)',&
                         & __LINE__,info)
                    GOTO 9999
                 ELSE
                    IF(ipackdata(1,1).NE.1.AND.ipackdata(1,1).NE.7.AND.    &
                         & ipackdata(1,1).NE.9.AND.ipackdata(1,1).NE.20) THEN
                       info = ppm_error_error
                       CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                            & 'ipackdata(1,1) must be element of {1,7,9,20}',&
                            & __LINE__,info)
                       GOTO 9999
                    END IF
                END IF
           END IF

           !--------------------------------------------------------------------
           ! check association of up, dup, bfr
           !--------------------------------------------------------------------
           IF(.NOT.ASSOCIATED(up)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'UP is empty',__LINE__,info)
              GOTO 9999
           END IF
           IF(.NOT.ASSOCIATED(dup)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'DUP is empty',__LINE__,info)
              GOTO 9999
           END IF
           IF(.NOT.ASSOCIATED(bfr)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'BFR is empty',__LINE__,info)
              GOTO 9999
           END IF
           
        END IF ! (ppm_debug.GT.0)

        mid = ppm_internal_mid(odeid)
        !-----------------------------------------------------------------------
        ! check state if finished, bail out
        !-----------------------------------------------------------------------
        IF(ppm_ode_state(mid).EQ.ppm_ode_state_finished) GOTO 9999
        !-----------------------------------------------------------------------
        ! check istage, if greater that ppm_ode_stages, then
        ! bail out
        !-----------------------------------------------------------------------
        IF(ppm_ode_stages(mid).LT.istage) GOTO 9999
        
        !-----------------------------------------------------------------------
        ! get times and scheme
        !-----------------------------------------------------------------------
        t      = time(3)
        dt     = time(4)
        IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
           scheme = ppm_ode_kscheme(mid)
        ELSE
           scheme = ppm_ode_ischeme(mid)
        END IF
        
        IF(ppm_ode_adaptive(mid)) THEN
           !--------------------------------------------------------------------
           ! compute adaptive timestep [TODO]
           !--------------------------------------------------------------------
           CALL ppm_error(ppm_err_argument,'ppm_ode_step',&
                & 'adaptivity not yet implemented',__LINE__,info)

           GOTO 9999
        END IF
        
        SELECT CASE(scheme)
        CASE(ppm_param_ode_scheme_eulerf)
           !--------------------------------------------------------------------
           !=======
           ! euler:
           !=======
           !--------------------------------------------------------------------

           !--------------------------------------------------------------------
           ! call right hand side and do an euler step
           !--------------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
#if      __MODE == __SCA           
           DO j=1,npart
              up(j) = up(j) + dt * dup(j)
           END DO
#elif    __MODE == __VEC
           DO i=1,lda
              DO j=1,npart
                 up(i,j) = up(i,j) + dt * dup(i,j)
              END DO
           END DO
#endif 
           t  = t  + dt

           IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
              ppm_ode_state(mid) = ppm_ode_state_running
           END IF
           ! how much to save of this
           ppm_ode_sent(mid) = 0

           CASE(ppm_param_ode_scheme_sts)
              !-----------------------------------------------------
              !  compute the new dt
              !-----------------------------------------------------
              stsn = ipackdata(1,1)
              if(present(rpackdata)) stsnu(stsn) = rpackdata(1,1)
              tau = dt/((stsnu(stsn)-1.0_mk)*&
             & COS((2.0_mk*REAL(istage,mk)-1.0_mk)/REAL(stsn,mk)*M_PI*0.5_mk)&
             & +1.0_mk+stsnu(stsn))
              !-----------------------------------------------------------------
              !=======
              !  euler:
              !=======
              !-----------------------------------------------------------------
              !  call right hand side and do an euler step
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
#if      __MODE == __SCA           
           DO j=1,npart
              up(j) = up(j) + tau * dup(j)
           END DO
#elif    __MODE == __VEC
           DO i=1,lda
              DO j=1,npart
                 up(i,j) = up(i,j) + tau * dup(i,j)
              END DO
           END DO
#endif
              t  = t  + tau
     
              IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
                 ppm_ode_state(mid) = ppm_ode_state_running
              END IF
              ! how much to save of this
              ppm_ode_sent(mid) = 0
           
        CASE(ppm_param_ode_scheme_tvdrk2)
           !--------------------------------------------------------------------
           !============
           ! 2nd tvd rk:
           !============
           !--------------------------------------------------------------------
           SELECT CASE(istage)
           CASE(1)
              !-----------------------------------------------------------------
              ! call rhs, save the old up and do an euler step
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
#if      __MODE == __SCA
              DO i=1,Npart
                 bfr(1,i) = up(i)
                 up(i) = up(i) + dt*dup(i)
              END DO
#elif    __MODE == __VEC
              DO i=1,Npart             
                 DO j=1,lda
                    bfr(j,i) = up(j,i)
                    up(j,i) = up(j,i) + dt*dup(j,i)
                 END DO
              END DO
#endif
              ppm_ode_sent(mid) = 1
           CASE(2) 
              !-----------------------------------------------------------------
              ! call rhs, and do another euler step
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
#if      __MODE == __SCA
              DO i=1,Npart
                 up(i) = up(i) + dt*dup(i)
              END DO
#elif    __MODE == __VEC
              DO i=1,Npart
                 DO j=1,lda
                    up(j,i) = up(j,i) + dt*dup(j,i)
                 END DO
              END DO
#endif
              !-----------------------------------------------------------------
              ! interpolate
              !-----------------------------------------------------------------
#if      __MODE == __SCA
              DO i=1,Npart
                 up(i)   = 0.5_mk * (up(i)   + bfr(1,i)    )
              END DO
#elif    __MODE == __VEC
              DO i=1,Npart                 
                 DO j=1,lda
                 up(j,i) = 0.5_mk * (up(j,i) + bfr(j,i))
              END DO
           END DO
#endif      
              t = t + dt
              ppm_ode_sent(mid) = 0
              IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
                 ppm_ode_state(mid) = ppm_ode_state_running
              END IF

           END SELECT
           
        CASE(ppm_param_ode_scheme_midrk2)
           !--------------------------------------------------------------------
           !=============
           ! mid point rk
           !=============
           !--------------------------------------------------------------------
           SELECT CASE(istage)
           CASE(1)
              !-----------------------------------------------------------------
              ! 
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"

#if      __MODE == __SCA
              DO i=1,Npart
                 bfr(1,i)     = up(i)
                 up(i) = up(i) + 0.5_mk* dt * dup(i)
              END DO
#elif    __MODE == __VEC
              DO j=1,lda
                 DO i=1,Npart           
                    bfr(j,i) = up(j,i)
                    up(j,i) = up(j,i) + 0.5_mk*dt*dup(j,i)
                 END DO
              END DO
#endif
              ppm_ode_sent(mid) = 1
           CASE(2) 
              !-----------------------------------------------------------------
              ! 
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
              !-----------------------------------------------------------------
              ! 
              !-----------------------------------------------------------------
#if      __MODE == __SCA
              DO i=1,Npart
                 up(i)   = bfr(1,i) + dt * dup(i)
              END DO
#elif    __MODE == __VEC
              DO i=1,Npart
                 DO j=1,lda
                    up(j,i) = bfr(j,i) + dt * dup(j,i)
                 END DO
              END DO
#endif                 
                 

              t = t + dt
              ppm_ode_sent(mid) = 0
              IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
                 ppm_ode_state(mid) = ppm_ode_state_running
              END IF

           END SELECT

           
        CASE(ppm_param_ode_scheme_rk4)
           !--------------------------------------------------------------------
           !=============
           ! Runge Kutta 4
           !=============
           !--------------------------------------------------------------------
           SELECT CASE(istage)
           CASE(1)
              !-----------------------------------------------------------------
              ! x_n + 1/2 dt k1
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
              DO i=1,Npart
#if      __MODE == __SCA
                 bfr(1,i)     =  up(i)
                 bfr(2,i)     = dup(i) ! k1
#elif    __MODE == __VEC             
                 bfr(1:lda,i) = up(:,i)
                 bfr((lda+1):2*lda,i) = dup(:,i) !k1
#endif
              END DO
              DO i=1,Npart
#if      __MODE == __SCA
                 up(i) = up(i) + 0.5_mk* dt * dup(i)
#elif    __MODE == __VEC
                 DO ilda=1,lda
                    up(ilda,i) = up(ilda,i) + 0.5_mk*dt*dup(ilda,i)
                 END DO
#endif       
              END DO
              ppm_ode_sent(mid) = 2
           CASE(2) 
              !-----------------------------------------------------------------
              ! 
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
              !-----------------------------------------------------------------
              ! x_n + 1/2 dt k2
              !-----------------------------------------------------------------
              DO i=1,Npart
#if      __MODE == __SCA
                 bfr(3,i)     = dup(i) !k2
#elif    __MODE == __VEC             
                 DO ilda=1,lda
                    bfr((2*lda+ilda),i) = dup(ilda,i) !k2
                 END DO
#endif                 

#if      __MODE == __SCA
                 up(i)   = bfr(1,i) + 0.5_mk* dt * dup(i)
#elif    __MODE == __VEC                 
                 DO ilda=1,lda
                    up(ilda,i) = bfr(ilda,i) + 0.5_mk*dt * dup(ilda,i)
                 END DO
#endif                 
                 
              END DO
              ppm_ode_sent(mid) = 3
           CASE(3) 
              !-----------------------------------------------------------------
              ! 
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
              !-----------------------------------------------------------------
              ! x_n + dt k3
              !-----------------------------------------------------------------
              DO i=1,Npart
#if      __MODE == __SCA
                 bfr(4,i)     = dup(i) !k3
#elif    __MODE == __VEC
                 DO ilda=1,lda
                    bfr((3*lda+ilda),i) = dup(ilda,i) !k3
                 END DO
#endif                 

#if      __MODE == __SCA
                 up(i)   = bfr(1,i) + dt * dup(i)
#elif    __MODE == __VEC                 
                 DO ilda=1,lda
                    up(ilda,i) = bfr(ilda,i) + dt * dup(ilda,i)
                 END DO
#endif                 
                 
              END DO
              ppm_ode_sent(mid) = 4
           CASE(4) 
              !-----------------------------------------------------------------
              ! 
              !-----------------------------------------------------------------
#include "ppm_ode_rhsfunc_macro.h"
              !-----------------------------------------------------------------
              ! x_n + 1/6 dt (k1 + 2 k2 + 2k3 +k4)
              !-----------------------------------------------------------------
              DO i=1,Npart
#if      __MODE == __SCA
                 up(i)   = bfr(1,i) + 1.0_MK/6.0_MK*dt*                        &
        &                  (bfr(2,i) + 2.0_MK*bfr(3,i) + 2.0_MK*bfr(4,i)+ dup(i))
#elif    __MODE == __VEC
        DO ilda=1,lda
           up(ilda,i) = bfr(ilda,i) +  1.0_MK/6.0_MK*dt*                   & 
           &                (bfr((lda+ilda),i) + 2.0_MK*bfr((2*lda+ilda),i)  &
           &                             + 2.0_MK*bfr((3*lda+ilda),i)+ dup(ilda,i))
        END DO
#endif                 

      END DO
              t = t + dt
              ppm_ode_sent(mid) = 0
              IF(ppm_ode_state(mid).EQ.ppm_ode_state_kickoff) THEN
                 ppm_ode_state(mid) = ppm_ode_state_running
              END IF

           END SELECT

        END SELECT
        !-----------------------------------------------------------------------
        ! pass back dt and time
        !-----------------------------------------------------------------------
        time(3) = t
        time(4) = dt
        
        !-----------------------------------------------------------------------
        ! stop ode if were ready
        !-----------------------------------------------------------------------
        IF(time(3).GE.time(2)) THEN
           ppm_ode_state(mid) = ppm_ode_state_finished
        END IF
        
9999    CONTINUE        
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop('ppm_ode_step',t0,info)
        RETURN
#if   __MODE == __SCA
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_ode_step_ss
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_ode_step_ds
#endif
#elif __MODE == __VEC
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_ode_step_sv
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_ode_step_dv
#endif
#endif
