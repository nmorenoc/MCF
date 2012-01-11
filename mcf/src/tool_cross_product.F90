      SUBROUTINE tool_cross_product(this,a,b,c,stat_info)
        !----------------------------------------------------
        ! Subroutine  : tool_corss_product
        !----------------------------------------------------
        !
        ! Purpose     : Calculate cross product of two
        !               input vectors
        !
        ! Routines    :
        !
        ! Remarks     : Input, output vector are 3D.
        !
        ! References  :
        !
        ! Revisions   : V0.1 05.10 2009, original version.
        !
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments     
        !----------------------------------------------------
        
        TYPE(Tool), INTENT(IN)                  :: this
        REAL(MK), DIMENSION(:),INTENT(IN)       :: a
        REAL(MK), DIMENSION(:),INTENT(IN)       :: b
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: c
        INTEGER, INTENT(OUT)                    :: stat_info

        INTEGER                                 :: num_dim
        
        
        !----------------------------------------------------
        ! Initialization
        !
        ! This is supposed to be used, otherwise,
        ! compiler complains that it is not used.
        !----------------------------------------------------
        
        stat_info = this%flag
        stat_info = 0
        
        
        num_dim = SIZE(a,1)
        IF ( num_dim /=3 ) THEN
           PRINT *, "tool_cross_product : ", &
                "a's dimesnion should be 3 !"
           stat_info = -1
           GOTO 9999
        END IF
        
        num_dim = SIZE(b,1)
        IF ( num_dim /=3 ) THEN
           PRINT *, "tool_cross_product : ", &
                "b's dimesnion should be 3 !"
           stat_info = -1
           GOTO 9999
        END IF
        
        num_dim = SIZE(c,1)
        IF ( num_dim /=3 ) THEN
           PRINT *, "tool_cross_product : ", &
                "c's dimesnion should be 3 !"
           stat_info = -1
           GOTO 9999
        END IF
        
        c(1) = a(2)*b(3)-a(3)*b(2)
        c(2) = a(3)*b(1)-a(1)*b(3)
        c(3) = a(1)*b(2)-a(2)*b(1)
        
        !----------------------------------------------------
        ! Return.
        !----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE tool_cross_product
      
      
