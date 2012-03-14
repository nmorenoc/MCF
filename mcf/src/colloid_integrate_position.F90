      SUBROUTINE colloid_integrate_position(this,dt,accuracy_order,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_integrate_position
        !----------------------------------------------------
        !
        ! Purpose     : Integrate the positions of centers of
        !               colloids.
        !
        ! Remark      : Colloid are modeled as rigid body.
        !
        !
        !  Revision   : V0.1  23.06.2009, original.
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
        
        TYPE(Colloid), INTENT(OUT)      :: this
        REAL(MK), INTENT(IN)            :: dt
        INTEGER, INTENT(IN)             :: accuracy_order
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

#ifdef __PARTICLES_POSITION_FIXED
#else
        IF( this%translate ) THEN
           
           !-------------------------------------------------
           ! Select different accuracy oder:
           ! 1 velocity contribution;
           ! 2 velocity + acceleration contribution.
           !-------------------------------------------------
           
           SELECT CASE (accuracy_order)
              
           CASE(1)
              
              this%x(:,:) = &
                   this%x(:, :) + this%v(:,:) * dt
              
           CASE(2)
           
              this%x(:,:) = &
                   this%x(:,:) + this%v(:,:) * dt + &
                   0.5_MK * this%f(:,:) * dt**2
              
           END SELECT ! accuracy_order
           
        END IF
     
#endif
   
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_integrate_position
      
