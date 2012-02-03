#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_map_push
      !-------------------------------------------------------------------------
      !
      !  Purpose      : pushes whats needed of the buffer
      !
      !  Input        : odeid                   (I) mode to push
      !                 bfr(:,:)               (F) buffer to push
      !                 lda                    (I) leading dimension
      !                 Npart                  (I) number of particles
      !
      !  Output       : info                   (I) return status
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_ode_map_push.f,v $
      !  Revision 1.8  2004/08/13 15:44:25  michaebe
      !  included mpart as input argument
      !
      !  Revision 1.7  2004/08/12 13:48:22  michaebe
      !  included check for ldasend -> bail out if 0
      !
      !  Revision 1.6  2004/08/12 13:10:39  michaebe
      !  corrected caller specification in substart
      !
      !  Revision 1.5  2004/07/26 13:49:19  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.4  2004/07/26 11:33:04  michaebe
      !  inserted the use of the ppm ode data module.
      !
      !  Revision 1.3  2004/07/26 07:50:51  michaebe
      !  Atomized. Otherwise no changes
      !
      !  Revision 1.2  2004/06/10 16:20:04  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.1  2004/02/19 08:33:56  michaebe
      !  initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_ode_map_push_s(odeid,bfr,lda,Npart,mpart,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_ode_map_push_d(odeid,bfr,lda,Npart,mpart,info)
#endif
        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
#include "ppm_define.h"

        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_map_part
        USE ppm_module_error
        USE ppm_module_substart
        USE ppm_module_data_ode
        USE ppm_module_substop
        USE ppm_module_data
        
        IMPLICIT NONE
#if     __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: mk = ppm_kind_single
#else
        INTEGER, PARAMETER :: mk = ppm_kind_double
#endif
        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER,                    INTENT(  out) :: info
        INTEGER,                    INTENT(in   ) :: odeid
        REAL(mk), DIMENSION(:,:),   POINTER       :: bfr
        INTEGER,                    INTENT(in   ) :: lda
        INTEGER,                    INTENT(in   ) :: Npart
        INTEGER,                    INTENT(inout) :: mpart

        !-----------------------------------------------------------------------
        ! Local Variables
        !-----------------------------------------------------------------------
        INTEGER                                   :: ldasend
        INTEGER                                   :: throwaway
        INTEGER                                   :: mid, umidmax, umidmin
        INTEGER                                   :: to_topo
        
        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart('ppm_ode_map_push',t0,info)
        
        !-----------------------------------------------------------------------
        ! general remark:
        ! were going to push only the stages that are needed to the
        ! other cpu. But that guy will pop the stuff and create a new
        ! array to comprise  only the stages that weve sent? so the
        ! array may shrink? must not happen. 
        !-----------------------------------------------------------------------
        
        IF(Npart.EQ.0) THEN
           !--------------------------------------------------------------------
           ! just save the number of stages that we would have sent
           !--------------------------------------------------------------------
           ! already happened in ppm_ode_step
           !--------------------------------------------------------------------
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        !  check input arguments
        !-----------------------------------------------------------------------
        IF(ppm_debug.GT.0) THEN
           !--------------------------------------------------------------------
           ! check if ppm is initialized
           !--------------------------------------------------------------------
           IF(.NOT.ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_ode_map_push',&
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
                   & 'odeid does not exist',__LINE__,info)
              GOTO 9999
           ELSE
              IF(ppm_internal_mid(odeid).EQ.-HUGE(odeid)) THEN
                 !--------------------------------------------------------------
                 ! user mid does not exist
                 !--------------------------------------------------------------
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_ode_step',& 
                      & 'odeid does not exist',__LINE__,info)
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
           IF(.NOT.ASSOCIATED(bfr)) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_step', &
                   & 'BFR is empty',__LINE__,info)
              GOTO 9999
           END IF
        END IF

        mid     = ppm_internal_mid(odeid)
        ldasend = lda*ppm_ode_sent(mid)
        IF(ldasend.EQ.0) GOTO 9999
        
        to_topo = -1
        CALL ppm_map_part(bfr,ldasend,Npart,mpart,to_topo,& 
             & ppm_param_map_push, info)
        IF(info.NE.0) THEN
           GOTO 9999
        END IF
        
        
        


9999    CONTINUE        
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop('ppm_ode_map_push',t0,info)
        RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_ode_map_push_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_ode_map_push_d
#endif
