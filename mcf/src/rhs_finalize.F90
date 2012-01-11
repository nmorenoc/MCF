      SUBROUTINE rhs_finalize(this,stat_info)
        
        TYPE(Rhs), INTENT(IN)                   :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = this%rhs_density_type
        stat_info = 0
        
        PRINT *, "rhs_finalize : ", "Finished!"
        
        RETURN
        
      END SUBROUTINE rhs_finalize
