!----------------------------------------------------
! Subroutine  : debug_get* routines.
!----------------------------------------------------
!
! Revisions   : V0.1 16.11.2010, original version.
!
!----------------------------------------------------
! Author      : Xin Bian
! Contact     : xin.bian@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!----------------------------------------------------      

      INTEGER FUNCTION debug_get_flag(this,stat_info)
        
        TYPE(Debug), INTENT(IN)         :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        debug_get_flag = this%flag
        
        RETURN
      END FUNCTION debug_get_flag
      

      REAL(MK) FUNCTION debug_get_time(this,stat_info)
        TYPE(Debug), INTENT(IN)         :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        REAL(MK)                        :: time

        stat_info = this%flag
        stat_info = 0
        
        CALL ppm_time(time,stat_info)
        

        debug_get_time = time
        
        RETURN
      END FUNCTION debug_get_time
