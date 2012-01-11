      SUBROUTINE control_finalize(this,stat_info)
        !----------------------------------------------------
        !  Subroutine   :  control_finalize
        !----------------------------------------------------
        !
        !  Purpose      : Destrutor of Class Control.
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : V0.1 01.03.2009, original version.
        !
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        TYPE(Control),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                   :: stat_info
        
        !------------------------
        ! Finalize
        !-----------------------
        
        stat_info = this%debug_flag
        stat_info = 0         

        PRINT *, "control_finalize : ", "Finished!"

        RETURN
        
      END SUBROUTINE control_finalize
     

