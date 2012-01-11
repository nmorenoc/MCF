      SUBROUTINE marching_finalize(this,stat_info)
          
        TYPE(Marching), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: stat_info_sub

        
        stat_info     = this%integrate_type
        stat_info     = 0
        stat_info_sub = 0
                
        CALL statistic_finalize(this%statis,stat_info_sub)
        
        PRINT *, "marching_finalize : ", "Finished!"

        RETURN
        
      END SUBROUTINE marching_finalize
