      SUBROUTINE statistic_finalize(this,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = this%num_dim
        stat_info = 0

        PRINT *, "statistic_finalize : ", "Finished!"
        
        RETURN          
        
      END SUBROUTINE statistic_finalize
      
