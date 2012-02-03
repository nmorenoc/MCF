#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_ode_map_pop
      !-------------------------------------------------------------------------
      !
      !  Purpose      : pops whats needed of the buffer
      !
      !  Input        : odeid                  (I) mode to push
      !                 lda                    (I) leading dimension
      !                 Npart                  (I) number of particles
      !
      !  Output       : bfr(:,:)               (F) buffer to pop into
      !                 info                   (I) RETURN status
      !
      !  Remarks      : 
      !
      !  References   : 
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_ode_map_pop.f,v $
      !  Revision 1.14  2006/09/04 18:34:53  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.13  2005/04/21 05:00:56  ivos
      !  Removed unnecessary USE ppm_module_error.
      !
      !  Revision 1.12  2004/08/13 15:44:25  michaebe
      !  included mpart as input argument
      !
      !  Revision 1.11  2004/08/13 12:27:06  michaebe
      !  corrected some mpart bug
      !
      !  Revision 1.10  2004/08/12 13:48:22  michaebe
      !  included check for ldasend -> bail out if 0
      !
      !  Revision 1.9  2004/08/12 13:10:58  michaebe
      !  corrected caller specificaion in substart
      !
      !  Revision 1.8  2004/07/26 13:49:18  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.7  2004/07/26 11:59:40  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.6  2004/07/26 11:47:17  michaebe
      !  forgot to use ppm_module_error
      !
      !  Revision 1.5  2004/07/26 11:33:04  michaebe
      !  inserted the use of the ppm ode data module.
      !
      !  Revision 1.4  2004/07/26 07:50:32  michaebe
      !  Atomized. Otherwise no changes.
      !
      !  Revision 1.3  2004/06/10 16:20:03  ivos
      !  Moved all cpp directtives to column 1. The NEC cpp did not recognize
      !  them otherwise!!!
      !
      !  Revision 1.2  2004/02/20 14:56:52  michaebe
      !  added the use of ppm_module_util.
      !
      !  Revision 1.1  2004/02/19 08:33:55  michaebe
      !  initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_ode_map_pop_s(odeid,bfr,lda,Npart,mpart,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_ode_map_pop_d(odeid,bfr,lda,npart,mpart,info)
#endif
        !-----------------------------------------------------------------------
        !  Includes
        !-----------------------------------------------------------------------
#include "ppm_define.h"

        !-----------------------------------------------------------------------
        !  Modules
        !-----------------------------------------------------------------------
        USE ppm_module_map_part
        USE ppm_module_alloc
        USE ppm_module_error
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_data_ode
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
        INTEGER,                    INTENT(inout) :: Mpart

        !-----------------------------------------------------------------------
        ! Local Variables
        !-----------------------------------------------------------------------
        INTEGER                                   :: ldasend
        INTEGER                                   :: throwaway
        INTEGER                                   :: mid, iopt
        INTEGER                                   :: to_topo, umidmin, umidmax
        INTEGER,                    DIMENSION(2)  :: dime

        !-----------------------------------------------------------------------
        !  call substart
        !-----------------------------------------------------------------------
        CALL substart('ppm_ode_map_pop',t0,info)
        
        IF(Mpart.EQ.0) THEN
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
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_ode_map_pop',&
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
              CALL ppm_error(ppm_err_argument,'ppm_ode_map_pop', &
                   & 'odeid does not exist',__LINE__,info)
              GOTO 9999
           ELSE
              IF(ppm_internal_mid(odeid).EQ.-HUGE(odeid)) THEN
                 !--------------------------------------------------------------
                 ! user mid does not exist
                 !--------------------------------------------------------------
                 info = ppm_error_error
                 CALL ppm_error(ppm_err_argument,'ppm_ode_map_pop',& 
                      & 'odeid does not exist',__LINE__,info)
                 GOTO 9999
              END IF
           END IF
           
           !--------------------------------------------------------------------
           ! check dimension
           !--------------------------------------------------------------------
           IF(Mpart.LT.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_map_pop', &
                   & 'Mpart cannot be <0',__LINE__,info)
              GOTO 9999
           END IF
           IF(lda.LE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_ode_map_pop', &
                   & 'LDA must be >00',__LINE__,info)
              GOTO 9999
           END IF
        END IF

        mid     = ppm_internal_mid(odeid)
        ldasend = lda*ppm_ode_sent(mid)
        IF(ldasend.EQ.0) GOTO 9999
        to_topo = -1
        !-----------------------------------------------------------------------
        ! get the stuff
        !-----------------------------------------------------------------------
        
        CALL ppm_map_part(bfr,ldasend,Npart,mpart,to_topo,& 
             & ppm_param_map_pop, info)
        IF(info.NE.0) THEN
           GOTO 9999
        END IF
        !-----------------------------------------------------------------------
        ! have to blow it up to full buffer size again 
        !-----------------------------------------------------------------------
        IF(ppm_ode_sent(mid).LT.ppm_ode_bfrsize(mid)) THEN
           iopt = ppm_param_alloc_fit_preserve
           dime(1) = ppm_ode_bfrsize(mid)*lda
           dime(2) = Mpart
           CALL ppm_alloc(bfr,dime,iopt,info)
           IF(info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_ode_map_pop', &
                   & 'growing buffer BFR',__LINE__,info)
              GOTO 9999
           END IF
           
        END IF
        

9999    CONTINUE        
        !-----------------------------------------------------------------------
        ! substop
        !-----------------------------------------------------------------------
        CALL substop('ppm_ode_map_pop',t0,info)
        RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_ode_map_pop_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_ode_map_pop_d
#endif
