      !-------------------------------------------------------------------------
      !  Subroutine   :                     debug_substop
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Called at the end of each subroutine. For timing
      !                 and debugging purposes.
      !
      !  Input        :  rank      debugging MPI rank
      !			 caller    (C) name of calling subroutine
      !
      !  Input/output : 
      !
      !  Output       : ctime        (F) current time
      !                 stat_info      (I) return status
      !
      !  Routines     : ppm_time
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      ! Author: contact
      !-------------------------------------------------------------------------
      !  
      !-------------------------------------------------------------------------

      SUBROUTINE debug_substop_p(this,rank,caller,time_start,stat_info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
        USE ppm_module_time    	  
      


      !-------------------------------------------------------------------------
      !  Includes 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
        TYPE(Debug), INTENT(IN)         :: this
        INTEGER, INTENT(IN)             :: rank
        CHARACTER(LEN=*), INTENT(IN)    :: caller
        REAL(MK), INTENT(IN)           :: time_start
        INTEGER,  INTENT(OUT)           :: stat_info
      !----------------------------------------------------------------------
      !  Local variables 
      !---------------------------------------------------------------------

        CHARACTER(LEN=MAX_CHAR)         :: cbuf
        REAL(MK)                        :: time_end
        INTEGER                         :: stat_info_sub

      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------

        stat_info     = 0
        stat_info_sub = 0
        time_end      =0.0

        CALL ppm_time(time_end,stat_info)

        WRITE(cbuf,'(A,E16.8,A)') 'Leaving & Took: ',time_end-time_start, " s"
        CALL debug_print_msg(this,rank,caller,cbuf,stat_info_sub)   
     

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------

9999    CONTINUE
      
        RETURN
      
      END SUBROUTINE debug_substop_p
    
      SUBROUTINE debug_substop_s(this,caller,time_start,stat_info)

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
        USE ppm_module_time    	  
      


      !-------------------------------------------------------------------------
      !  Includes 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
        TYPE(Debug), INTENT(IN)         :: this
        CHARACTER(LEN=*), INTENT(IN)    :: caller
        REAL(MK), INTENT(IN)            :: time_start
        INTEGER,  INTENT(OUT)           :: stat_info
      !----------------------------------------------------------------------
      !  Local variables 
      !---------------------------------------------------------------------

        INTEGER                         :: rank
        CHARACTER(LEN=MAX_CHAR)         :: cbuf
        REAL(MK)                        :: time_end
        INTEGER                         :: stat_info_sub

      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------

        rank          = -1
        stat_info     = 0
        stat_info_sub = 0
        time_end      = 0.0

        CALL ppm_time(time_end,stat_info)
        
        WRITE(cbuf,'(A,E16.8,A)') 'Leaving & Took: ',time_end-time_start, " s"
        CALL debug_print_msg(this,rank,caller,cbuf,stat_info_sub)   
     

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------

9999    CONTINUE
      
        RETURN
      
      END SUBROUTINE debug_substop_s
    
