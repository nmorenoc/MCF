      SUBROUTINE stateEquation_finalize(this,stat_info)
  
        TYPE(StateEquation), INTENT(IN)         :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = this%stateEquation_type
        stat_info = 0
        
        PRINT *, "stateEquation_finalize : ", "Finished!"
        
        RETURN
        
      END SUBROUTINE stateEquation_finalize
