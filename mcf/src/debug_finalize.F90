      SUBROUTINE debug_finalize(this,stat_info)

        TYPE(Debug), INTENT(IN)         :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = this%flag
        stat_info = 0          
        
        RETURN

      END SUBROUTINE debug_finalize
