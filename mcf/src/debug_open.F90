      SUBROUTINE debug_open_time(this,rank,stat_info)
        !----------------------------------------------------
        ! Subroutine  : debug_open_time
        !----------------------------------------------------
        !
        ! Purpose     : Open time file for writting.
        !
        ! Input       : rank  : MPI rank
        !
        ! Output      : stat_info : return status
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
        INTEGER, INTENT(IN)             :: rank
        INTEGER,  INTENT(OUT)           :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        IF ( rank /=0 ) THEN
           PRINT *, "debug_open_time : ", &
                "rank must be 0 !"
           stat_info = -1
           GOTO 9999
        END IF
        
        OPEN(UNIT=this%time_file_unit,FILE=this%time_file,&
             STATUS="NEW",FORM="FORMATTED",ACTION="WRITE", &
             IOSTAT=stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "debug_open_time : ", &
                "Opening file ", &
                this%time_file, " has problem !"
           stat_info = -1
           GOTO 9999
        ELSE
           PRINT *,"New debug time file ", &
                TRIM(this%time_file), " being created !"
        END IF

9999    CONTINUE
        
        RETURN
      END SUBROUTINE debug_open_time
