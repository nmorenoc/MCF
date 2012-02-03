        !-----------------------------------------------------------------------
        ! ppm_ode_ischeme
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_ischeme, tlda, iopt, info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc, 'ppm_ode_modalloc.h', &
                &'PPM_ODE_ISCHEME',__LINE__,info)
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! ppm_ode_adaptive
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_adaptive, tlda, iopt, info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc, 'ppm_ode_modalloc.h', &
                &'PPM_ODE_ADAPTIVE',__LINE__,info)
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! ppm_ode_stages
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_stages, tlda, iopt, info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc, 'ppm_ode_modalloc.h', &
                &'PPM_ODE_STAGES',__LINE__,info)
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! ppm_ode_state
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_state, tlda, iopt, info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc, 'ppm_ode_modalloc.h', &
                &'PPM_ODE_STATE',__LINE__,info)
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! ppm_ode_sent
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_sent, tlda, iopt, info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc, 'ppm_ode_modalloc.h', &
                &'PPM_ODE_SENT',__LINE__,info)
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! ppm_ode_bfrsize
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_bfrsize, tlda, iopt, info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc, 'ppm_ode_modalloc.h', &
                &'PPM_ODE_BFRSIZE',__LINE__,info)
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! ppm_ode_sent
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_sent, tlda, iopt, info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc, 'ppm_ode_modalloc.h', &
                &'PPM_ODE_SENT',__LINE__,info)
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! ppm_internal_mid
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_internal_mid, tlda, iopt, info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc, 'ppm_ode_modalloc.h', &
                &'PPM_INTERNAL_MID',__LINE__,info)
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! ppm_user_mid
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_user_mid, tlda, iopt, info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc, 'ppm_ode_modalloc.h', &
                &'PPM_USER_MID',__LINE__,info)
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! ppm_ode_kscheme
        !-----------------------------------------------------------------------
        CALL ppm_alloc(ppm_ode_kscheme, tlda, iopt, info)
        IF(info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc, 'ppm_ode_modalloc.h', &
                &'PPM_ODE_KSCHEME',__LINE__,info)
           GOTO 9999
        END IF
