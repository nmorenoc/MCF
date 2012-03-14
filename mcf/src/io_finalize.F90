      SUBROUTINE io_finalize(this,rank,num_colloid,stat_info)

        TYPE(IO),INTENT(IN)             :: this
        INTEGER, INTENT(IN)             :: rank
        INTEGER, INTENT(IN)             :: num_colloid
        INTEGER,INTENT(OUT)             :: stat_info
        
        INTEGER                         :: i
        LOGICAL                         :: lopened
        
        stat_info = 0
        
        IF ( rank /= 0 ) THEN
           PRINT *, __FILE__, ":", __LINE__
           stat_info = -1
           GOTO 9999
        END IF
        
        INQUIRE(UNIT=this%statistic_unit,OPENED=lopened)
        
        IF ( lopened ) THEN  
           CLOSE(this%statistic_unit)
        END IF
        
        DO i = 1, num_colloid
           
           INQUIRE(UNIT=this%colloid_unit+i,OPENED=lopened)
           
           IF ( lopened ) THEN           
              CLOSE(this%colloid_unit)
           END IF
           
        END DO
        
        INQUIRE(UNIT=this%boundary_unit,OPENED=lopened)
        
        IF ( lopened ) THEN  
           CLOSE(this%boundary_unit)
        END IF
        
        PRINT *, "io_finalize : ", "Finished!"
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE io_finalize
