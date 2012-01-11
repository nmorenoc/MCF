      SUBROUTINE debug_init(this,d_flag,stat_info)

        TYPE(Debug), INTENT(OUT)        :: this
        INTEGER, INTENT(IN)             :: d_flag
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%flag      = d_flag
        this%time_file = TRIM("debug_time.dat")
        this%time_file_unit = 1000
        
        RETURN
        
      END SUBROUTINE debug_init
      

      SUBROUTINE debug_display_parameters(this,stat_info)
        
        TYPE(Debug),INTENT(IN)          ::this
        INTEGER,INTENT(OUT)             :: stat_info
        
        stat_info = 0
        
        PRINT *, '------------------Start------------------'
        PRINT *, '     Debug parameters'
        PRINT *, '-----------------------------------------'
        
                
        PRINT *, "Debug flag : ", this%flag
        
        PRINT *, '-------------------End-------------------'

        RETURN          
        
      END SUBROUTINE debug_display_parameters
      
