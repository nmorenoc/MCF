      SUBROUTINE kernel_finalize(this,stat_info)
        
        TYPE(Kernel), INTENT(IN)        :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = this%kernel_type
        stat_info = 0

        PRINT *, "kernel_finalize : ", "Finished!"
        
        RETURN
        
      END SUBROUTINE kernel_finalize
      
