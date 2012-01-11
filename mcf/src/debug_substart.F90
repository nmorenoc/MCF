      SUBROUTINE debug_substart_p(this,rank,caller,time_start,stat_info)
        !----------------------------------------------------
        ! Subroutine  : debug_substart
        !----------------------------------------------------
        !
        ! Purpose     : Called at the start of calling
        !               subroutine showing entering message.
        !
        ! Input       : rank  : MPI rank
        !		caller: name of calling subroutine
        !
        ! Output      : ctime     : current time
        !               stat_info : return status
        !
        ! Routines    : ppm_time
        !
        ! Remarks     :
        !
        ! References  :
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
        
        !----------------------------------------------------
        ! Arguments 
        !----------------------------------------------------
        
        TYPE(Debug), INTENT(IN)         :: this
        INTEGER, INTENT(IN  )           :: rank
        CHARACTER(LEN=*), INTENT(IN)    :: caller
        REAL(MK), INTENT(OUT)           :: time_start
        INTEGER,  INTENT(OUT)           :: stat_info
        
        !----------------------------------------------------
        ! Local variables 
        !----------------------------------------------------
        
        INTEGER                           :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        time_start    =0.0
        
        CALL ppm_time(time_start,stat_info)
        
        CALL debug_print_msg(this,rank,caller,'Entering.',stat_info_sub)
        
        !----------------------------------------------------
        ! Return 
        !----------------------------------------------------

9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE debug_substart_p
      
      
      SUBROUTINE debug_substart_s(this,caller,time_start,stat_info)
        
        !----------------------------------------------------
        ! Arguments 
        !----------------------------------------------------
        
        TYPE(Debug), INTENT(IN)         :: this
        CHARACTER(LEN=*), INTENT(IN)    :: caller
        REAL(MK), INTENT(OUT)           :: time_start
        INTEGER,  INTENT(OUT)           :: stat_info
        
        !----------------------------------------------------
        ! Local variables 
        !----------------------------------------------------
        
        INTEGER                         :: rank
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------
        
        rank          = -1
        stat_info     = 0
        stat_info_sub = 0
        time_start    = 0.0
        
        CALL ppm_time(time_start,stat_info)
        
        CALL debug_print_msg(this,rank,caller,'Entering.',stat_info_sub)
        
        
        !----------------------------------------------------
        ! Return 
        !----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE debug_substart_s
      
