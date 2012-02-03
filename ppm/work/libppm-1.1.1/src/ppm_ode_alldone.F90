#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Function     :                 ppm_ode_alldone
      !-------------------------------------------------------------------------
      !
      !  Purpose      : will say .true. if all modes are integrated 
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
      !  $Log: ppm_ode_alldone.f,v $
      !  Revision 1.6  2004/07/26 11:59:40  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.5  2004/07/26 11:33:03  michaebe
      !  inserted the use of the ppm ode data module.
      !
      !  Revision 1.4  2004/07/26 11:28:13  michaebe
      !  syntax error... removed.
      !
      !  Revision 1.3  2004/07/26 07:46:51  michaebe
      !  Atomized. Otherwise no changes.
      !
      !  Revision 1.2  2004/06/10 16:20:03  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.1  2004/02/19 08:33:53  michaebe
      !  initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      FUNCTION ppm_ode_alldone(info)
        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
#include "ppm_define.h"

        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_error
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_data_ode
        USE ppm_module_data
        IMPLICIT NONE

        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER,                    INTENT(  OUT) :: info
        LOGICAL                                   :: ppm_ode_alldone
        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart('ppm_ode_alldone',t0,info)

        IF(ppm_debug.GT.0) THEN
           IF(.NOT.ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_ode_alldone',&
                   & 'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
           END IF
        END IF
        
        IF(ANY(ppm_ode_state(1:ppm_max_mid) .NE. ppm_ode_state_finished)) THEN
           ppm_ode_alldone = .FALSE.
        ELSE
           ppm_ode_alldone = .TRUE.
        END IF
        
        

9999    CONTINUE        
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop('ppm_ode_alldone',t0,info)
        RETURN
      END FUNCTION ppm_ode_alldone
