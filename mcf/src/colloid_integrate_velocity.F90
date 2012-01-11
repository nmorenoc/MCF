      SUBROUTINE colloid_integrate_velocity(this,dt,lamda,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_integrate_velocity
        !----------------------------------------------------
        !
        ! Purpose     : Integrate translation velocity of 
        !               centers of colloid, also
        !               its rotation velocity.
        !
        ! Remark      : Colloid are modelled as rigid body.
        !               
        !
        ! Revision    : V0.2 05.10.2009, including rotating
        !               velocity.
        !
        !               V0.1 23.06.2009,original version.
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
        ! Arguments:
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(OUT)      :: this
        REAL(MK), INTENT(IN)            :: dt
        REAL(MK), INTENT(IN)            :: lamda
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables:
        !----------------------------------------------------

        INTEGER                         :: dim
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info  = 0
        dim        = this%num_dim
        
        
        IF ( this%translate ) THEN
            
           this%v(:,:) = &
                this%v(:,:) + &
                lamda * this%f(:,:) * dt
           
        END IF ! translate
        
        IF ( this%rotate ) THEN
           
           this%omega(:,:) = &
                this%omega(:,:) + &
                lamda * this%alpha(:,:) * dt
           
        END IF ! rotate
        
        RETURN
        
      END SUBROUTINE colloid_integrate_velocity
      
