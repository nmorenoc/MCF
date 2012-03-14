      SUBROUTINE colloid_in_relative_position(this,x,sid, &
           rx,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_in_relative_position
        !----------------------------------------------------
        ! Purpose     : Compute boundary particle's  
        !               position relative to the center
        !               of colloid, according to different 
        !               boundary conditions.
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     : x the absolute position of a boundary
        !               particle in the simulation. sid is its
        !               species ID. rx is the relative position
        !               of boundary particle to colloid center,
        !               taking into account images of colloid
        !               by boundary condition.
        !
        ! Revisions   : V0.1 23.08.2010, original version
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
        !
        ! Input 
        !
        ! this  : object of a colloid.
        ! x     : position of a fluid particle.
        ! sid   : species ID of a colloid boundary particle.
        !
        ! Output
        !
        ! rx    : relative position of the boundary particle
        !         to the center of the colloid.
        !         
        ! stat_info : status of the routine.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: x
        INTEGER, INTENT(IN)                     :: sid
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: rx
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local parameters:
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        REAL(MK), DIMENSION(3)                  :: xcoll
        REAL(MK), DIMENSION(3)                  :: vcoll

        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        dim           = this%num_dim
        
        CALL colloid_in_nearest_image(this, x(1:dim), sid, &
             xcoll(1:dim),rx(1:dim),vcoll(1:dim),stat_info_sub)
        
        IF ( stat_info_sub /=0 ) THEN
           
           PRINT *, "colloid_in_relative_position: ", &
                "colloid_in_nearst_image failed"
           stat_info = -1
           GOTO 9999
           
        END IF
        
9999    CONTINUE
        
        
        RETURN
        
      END SUBROUTINE colloid_in_relative_position
      
