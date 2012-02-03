#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_create_ode
      !-------------------------------------------------------------------------
      !
      !  Purpose      : creates a mode using a given schemes
      !
      !  Input        : ischeme                (I) integration scheme
      !                 kscheme                (I) kickoff scheme
      !                 adaptive               (L) adaptive dt flag
      !
      !  Input/Output : odeid                  (I) user mode/ode id      
      !
      !  Output       : info                   (I) return status
      !                 bfrsize                (I) size of the buffer the 
      !                                            schemes need (user must
      !                                            allocate lda*bfrsize
      !                 nstage                 (I) number of stages
      ! 
      !  Remarks      : creates an odeid if odeid=-1
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_ode_create_ode.f,v $
      !  Revision 1.9  2006/09/04 18:34:53  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.8  2005/08/05 09:27:32  ivos
      !  Added STS stuff (credits: MB)
      !
      !  Revision 1.7  2004/08/12 13:24:08  michaebe
      !  removed print statement
      !
      !  Revision 1.6  2004/07/26 15:38:51  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.5  2004/07/26 11:33:03  michaebe
      !  inserted the use of the ppm ode data module.
      !
      !  Revision 1.4  2004/07/26 07:47:31  michaebe
      !  Atomized. Major bug-fixes.
      !
      !  Revision 1.3  2004/06/10 16:20:03  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.2  2004/02/20 14:55:38  michaebe
      !  added the use of ppm_module_util.
      !
      !  Revision 1.1  2004/02/19 08:33:54  michaebe
      !  initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_ode_create_ode(odeid,bfrsize,nstage,ischeme,kscheme,&
           &                        adaptive,info)
        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
#include "ppm_define.h"
      
        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_util_invert_list
        USE ppm_module_error
        USE ppm_module_alloc
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_data_ode
        USE ppm_module_data
                
        IMPLICIT NONE
        
        !-----------------------------------------------------------------------
        !  Arguments
        !-----------------------------------------------------------------------
        INTEGER,                    INTENT(  out) :: info
        INTEGER,                    INTENT(inout) :: odeid
        INTEGER,                    INTENT(  out) :: bfrsize
        INTEGER,                    INTENT(  out) :: nstage
        INTEGER,                    INTENT(in   ) :: ischeme
        INTEGER,  OPTIONAL,         INTENT(in   ) :: kscheme
        LOGICAL,  OPTIONAL,         INTENT(in   ) :: adaptive

        INTEGER                                   :: mid, tkscheme, iopt
        INTEGER, DIMENSION(3)                     :: tlda
        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart('ppm_ode_create_ode',t0,info)
        CALL ppm_module_data_ode_activate
        !-----------------------------------------------------------------------
        !  check input arguments
        !-----------------------------------------------------------------------
        IF(ppm_debug.GT.0) THEN
           !--------------------------------------------------------------------
           ! check if ppm is initialized
           !--------------------------------------------------------------------
           IF(.NOT.ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_ode_create_ode',&
              & 'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
           END IF
           !--------------------------------------------------------------------
           ! check that schemes are valid
           !--------------------------------------------------------------------
           IF(ischeme.GT.0.AND.ischeme.LE.SIZE(ppm_ode_scheme_o)) THEN
              IF(ppm_ode_scheme_o(ischeme).LT.1) THEN
                 !--------------------------------------------------------------
                 ! scheme not yet implemented
                 !--------------------------------------------------------------
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_ode_create_ode',&
                      & 'ISCHEME not implemented',__LINE__,info)
                 GOTO 9999
              END IF
           ELSE
              !-----------------------------------------------------------------
              ! scheme does not exist
              !-----------------------------------------------------------------
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_create_ode',&
                   & 'ISCHEME not implemented',__LINE__,info)
              GOTO 9999
           END IF
           IF(PRESENT(kscheme)) THEN
              IF(kscheme.GT.0.AND.kscheme.LE.SIZE(ppm_ode_scheme_o)) THEN
                 IF(ppm_ode_scheme_o(kscheme).LT.1) THEN
                    !-----------------------------------------------------------
                    ! scheme not yet implemented
                    !-----------------------------------------------------------
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_ode_create_ode',&
                         & 'KSCHEME not implemented',__LINE__,info)
                    GOTO 9999
                 END IF
              ELSE
                 !--------------------------------------------------------------
                 ! scheme does not exist
                 !--------------------------------------------------------------
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_ode_create_ode',&
                      & 'KSCHEME not implemented',__LINE__,info)
                 GOTO 9999
              END IF
           END IF
           !--------------------------------------------------------------------
           ! check odeid
           !--------------------------------------------------------------------
           IF(odeid.GE.0) THEN
              IF(ANY(ppm_user_mid(1:ppm_max_mid).EQ.odeid)) THEN
                 !--------------------------------------------------------------
                 ! already exists. sorry
                 !--------------------------------------------------------------
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_ode_create_ode',&
                      & 'ODEID already exists',__LINE__,info)
                 GOTO 9999
              END IF
           END IF
        END IF
        
        !-----------------------------------------------------------------------
        ! create or use a kickoff scheme
        !-----------------------------------------------------------------------
        IF(PRESENT(kscheme)) THEN
           tkscheme = kscheme
        ELSE
           tkscheme = ppm_ode_scheme_k(ischeme)
        END IF
        
        !-----------------------------------------------------------------------
        ! check if there is enough space in the arrays
        !-----------------------------------------------------------------------
        mid = ppm_max_mid + 1
        IF(mid.GT.ppm_max_mid_allocd) THEN
           ppm_max_mid_allocd = ppm_max_mid_allocd + 1
           iopt = ppm_param_alloc_grow_preserve
           tlda(1) = ppm_max_mid_allocd
#include "ppm_ode_modalloc.h"
        END IF
        
        !-----------------------------------------------------------------------
        ! if odeid.lt.0 then generate it for the user
        !-----------------------------------------------------------------------
        IF(odeid.LT.0) THEN
           odeid = mid
        END IF
        ppm_user_mid(mid) = odeid

        !-----------------------------------------------------------------------
        ! update the inverse list
        !-----------------------------------------------------------------------
        CALL ppm_util_invert_list(ppm_user_mid,ppm_internal_mid,info)
        ppm_max_mid = maxval(ppm_internal_mid)

        !-----------------------------------------------------------------------
        ! store stuff
        !-----------------------------------------------------------------------
        ppm_ode_ischeme(mid)     = ischeme
        ppm_ode_kscheme(mid)     = tkscheme
        IF(PRESENT(adaptive)) THEN
           ppm_ode_adaptive(mid) = adaptive
        ELSE
           ppm_ode_adaptive(mid) = .FALSE.
        END IF
        ppm_ode_state(mid)       = ppm_ode_state_inited
        ppm_ode_stages(mid)      = MAX(ppm_ode_scheme_s(ischeme),&
             &                         ppm_ode_scheme_s(kscheme))
        ppm_ode_bfrsize(mid)     = MAX(ppm_ode_scheme_m(ischeme),&
             &                         ppm_ode_scheme_m(kscheme))
        bfrsize = ppm_ode_bfrsize(mid)
        nstage  = ppm_ode_stages(mid)
9999    CONTINUE        
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop('ppm_ode_create_ode',t0,info)
        RETURN
      END SUBROUTINE ppm_ode_create_ode
