      SUBROUTINE tool_rotation_vector(this,dim,rm,axis,phi,stat_info)
        !----------------------------------------------------
        ! Subroutine  : tool_rotation_vector
        !----------------------------------------------------
        !
        ! Purpose     : Calculate rotation vector, i.e.,
        !               rotating axis unit vector I and 
        !               angle about I, given a rotation matrix.
        !
        ! Routines    :
        !
        ! Remarks     : dim: specify problem dimension;
        !               rm : rotation matrix.
        !               ax : Rotating axis is always 3D;
        !               phi: rotate angle;
        !
        ! References  : Chen et al. 
        !               Physics of Fluids 18, 103605, 2006
        !
        ! Revisions   : V0.1 23.08 2010, original version.
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
        INTEGER,INTENT(IN)                      :: dim
        REAL(MK),DIMENSION(:,:),INTENT(IN)      :: rm
        REAL(MK),DIMENSION(:),INTENT(OUT)       :: axis
        REAL(MK),INTENT(OUT)                    :: phi
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim1
        
        !----------------------------------------------------
        ! Initialization
        !
        ! This is supposed to be used, otherwise,
        ! compiler complains that it is not used.
        !----------------------------------------------------
        
        stat_info = this%flag
        stat_info = 0
        
        IF ( dim /=2 .AND. dim /=3 ) THEN
           PRINT *, "tool_rotation_vector: ", &
                "Dimension should be either 2 or 3!"
           stat_info = -1
           GOTO 9999
        END IF
        
        dim1 = SIZE(rm,1)
        IF ( dim1 /= 3 ) THEN
           PRINT *, "tool_rotation_vector: ", &
                "dimension of rotation matrix does not match!"
           stat_info = -1
           GOTO 9999
        END IF
        
        dim1 = SIZE(axis,1)
        IF ( dim1 /= 3 ) THEN
           PRINT *, "tool_rotation_vector: ", &
                "rotation axis should be 3D!"
           stat_info = -1
           GOTO 9999
        END IF
       
        
        phi = ACOS((rm(1,1)+rm(2,2)+rm(3,3)-1)/2.0_MK)
        axis(1) = ( rm(3,2) - rm(2,3) ) / (2.0_MK * SIN(phi))
        axis(2) = ( rm(1,3) - rm(3,1) ) / (2.0_MK * SIN(phi))
        axis(3) = ( rm(2,1) - rm(1,2) ) / (2.0_MK * SIN(phi))
        
        !----------------------------------------------------
        ! Return.
        !----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE tool_rotation_vector
      
      
