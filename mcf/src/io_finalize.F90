      SUBROUTINE io_finalize(this,num_colloid,stat_info)

        TYPE(IO),INTENT(IN)             :: this
        INTEGER, INTENT(IN)             :: num_colloid
        INTEGER,INTENT(OUT)             :: stat_info
        
        INTEGER                         :: i
        LOGICAL                         :: lopened
        
        stat_info = 0
        
        INQUIRE(UNIT=this%statistic_unit,OPENED=lopened)
        
        IF (lopened ) THEN  
           CLOSE(this%statistic_unit)
        END IF
        
        INQUIRE(UNIT=this%colloid_unit,OPENED=lopened)
        IF (lopened ) THEN           
           CLOSE(this%colloid_unit)
        END IF
        
        PRINT *, "io_finalize : ", "Finished!"

        RETURN          
        
      END SUBROUTINE io_finalize
