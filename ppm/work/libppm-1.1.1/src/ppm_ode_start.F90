#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_start
      !-------------------------------------------------------------------------
      !
      !  Purpose      : (re)starts the ode solver
      !
      !  Input        : 
      !
      !  Output       : info                   (I) return status
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_ode_start.f,v $
      !  Revision 1.8  2005/08/05 09:24:49  ivos
      !  restored corrected subroutine name.
      !
      !  Revision 1.7  2005/08/05 09:23:49  ivos
      !  Restored STS stuff from dom.
      !
      !  Revision 1.5  2004/07/26 13:49:19  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.4  2004/07/26 11:33:04  michaebe
      !  inserted the use of the ppm ode data module.
      !
      !  Revision 1.3  2004/07/26 07:53:25  michaebe
      !  Atomized. Inserted ppm_init check. Otherwise no changes.
      !
      !  Revision 1.2  2004/06/10 16:20:04  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
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
      SUBROUTINE ppm_ode_start(info)
        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
#include "ppm_define.h"

        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_error
        USE ppm_module_data_ode
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_data
                
        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER,                    INTENT(  out) :: info
        
        !-----------------------------------------------------------------------
        ! Local Variables
        !-----------------------------------------------------------------------
        INTEGER                                   :: mid
        
        !-----------------------------------------------------------------------
        ! call substart
        !-----------------------------------------------------------------------
        CALL substart('ppm_ode_start',t0,info)

        !-----------------------------------------------------------------------
        ! check stuff
        !-----------------------------------------------------------------------
        IF(ppm_debug.GT.0) THEN
           IF(.NOT.ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_ode_start',&
                   & 'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
           END IF
        END IF
        
        !-----------------------------------------------------------------------
        ! check if everybody is ready to run
        !-----------------------------------------------------------------------
        DO mid=1,ppm_max_mid
           IF(ppm_ode_state(mid).NE.ppm_ode_state_inited) THEN
              IF(ppm_ode_state(mid).NE.ppm_ode_state_finished) THEN
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_notready,'ppm_ode_start', &
                      'some ODEs not ready',__LINE__,info)
                 GOTO 9999
              END IF
           END IF
        END DO
        
        !-----------------------------------------------------------------------
        ! yes, so set them into kickoff state
        !-----------------------------------------------------------------------
        DO mid=1,ppm_max_mid
           ppm_ode_state(mid) = ppm_ode_state_kickoff
        END DO
        
9999    CONTINUE        
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop('ppm_ode_start',t0,info)
        RETURN
      END SUBROUTINE ppm_ode_start
