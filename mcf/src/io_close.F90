!------------------------------------------------------------
!  SUBROUTINE : close files to finish writing
!------------------------------------------------------------
      SUBROUTINE io_close(this,&
           num_statis,num_boundary,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_close
        !----------------------------------------------------
        !
        ! Purpose     : Close files for finishing writting.
        !
        ! Reference   :
        !
        ! Remark      : 
        !
        ! Revisions   : V0.1 07.12 2009, original version.
        !
        !----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments :
        !----------------------------------------------------
        
        TYPE(IO), INTENT(IN)            :: this
        INTEGER,INTENT(IN)              :: num_statis
        INTEGER,INTENT(IN)              :: num_boundary
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_colloid

        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
                
        IF ( num_statis > 0 ) THEN
           
           CALL io_close_statistic(this,stat_info_sub)
           
        END IF

        IF ( num_boundary > 0 ) THEN

           CALL io_close_boundary(this,stat_info_sub)
           
        END IF

#ifdef __IO_COLLOID_SEPARATE

        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF ( num_colloid > 0 ) THEN
           
           CALL io_close_colloid(this,num_colloid,stat_info_sub)
           
        END IF
        
#endif
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE io_close

      
      SUBROUTINE io_close_statistic(this,stat_info)
        TYPE(IO), INTENT(IN)         :: this
        INTEGER, INTENT(OUT)            :: stat_info

        LOGICAL                         :: lopened
        
        stat_info = 0
        
        INQUIRE(UNIT=this%statistic_unit,OPENED=lopened)
        
        IF (lopened ) THEN           
           CLOSE(this%statistic_unit)
        END IF
        
      END SUBROUTINE io_close_statistic

      
      SUBROUTINE io_close_statistic_relax(this,stat_info)
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        LOGICAL                         :: lopened
        
        stat_info = 0
        
        INQUIRE(UNIT=this%statistic_relax_unit,OPENED=lopened)
        
        IF (lopened ) THEN           
           CLOSE(this%statistic_relax_unit)
        END IF
        
      END SUBROUTINE io_close_statistic_relax
      
      
      SUBROUTINE io_close_colloid(this,num_colloid,stat_info)
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)             :: num_colloid
        INTEGER, INTENT(OUT)            :: stat_info

        INTEGER                         :: i
        LOGICAL                         :: lopened


        stat_info = 0
        
        DO i = 1, num_colloid
           
           INQUIRE(UNIT=this%colloid_unit+i,OPENED=lopened)
           IF (lopened ) THEN           
              CLOSE(this%colloid_unit+i)
           END IF
        END DO
        
      END SUBROUTINE io_close_colloid
      
      
      SUBROUTINE io_close_boundary(this,stat_info)
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(OUT)            :: stat_info

        LOGICAL                         :: lopened
        
        stat_info = 0
        
        INQUIRE(UNIT=this%boundary_unit,OPENED=lopened)
        
        IF (lopened ) THEN           
           CLOSE(this%boundary_unit)
        END IF
        
      END SUBROUTINE io_close_boundary

