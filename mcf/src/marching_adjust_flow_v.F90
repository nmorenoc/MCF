      SUBROUTINE marching_adjust_flow_v(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : marching_adjust_flow_v
        !----------------------------------------------------
        !
        ! Purpose     : Adjusting flow velocity.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 07.12 2009, original version.
        !
        !               V0.2 28.06.2010, pointed out by
        !               Adolfo that change of df is not
        !               a good idea.
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
        ! this       : an object of Marching Class.
        ! stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Marching), INTENT(INOUT)   :: this
        INTEGER, INTENT(OUT)            :: stat_info
        

	!----------------------------------------------------
      	! Local variables starts here :
      	!----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Control parameters :
        !
     	!----------------------------------------------------
      
        LOGICAL                         :: flow_v_fixed

        !----------------------------------------------------
        ! Physics parameters :
        !
        ! below : indicate if curret flow velocity is below 
        !         desired flow velocity, in case we need 
        !         constant flow velocity at far field.
     	!----------------------------------------------------
    
        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(:), POINTER :: body_force
        REAL(MK), DIMENSION(:), POINTER :: body_force_d
        INTEGER                         :: fd
        REAL(MK)                        :: flow_v                
        REAL(MK), DIMENSION(:), POINTER :: flow_v_current

	!----------------------------------------------------
        ! Colloid parameters :
     	!----------------------------------------------------
        
        INTEGER                         :: num_colloid
        TYPE(Colloid), POINTER          :: colloids
                
        

        !----------------------------------------------------
        ! Initialization of variables.
      	!----------------------------------------------------
       
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(body_force)
        NULLIFY(body_force_d)
        NULLIFY(flow_v_current)
        
        num_colloid = 0
        NULLIFY(colloids)
     
        !----------------------------------------------------
        ! Get control parameters.
        !----------------------------------------------------
        
        flow_v_fixed = &
             control_get_flow_v_fixed(this%ctrl,stat_info_sub)
        
        IF ( .NOT. flow_v_fixed ) THEN
           PRINT *, "marching_adjust_flow_v : ", &
                "flow velocity is not fixed ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Get physics parameters.
        !----------------------------------------------------
        
        num_dim = &
             physics_get_num_dim(this%phys,stat_info_sub)
        
        flow_v = &
             physics_get_flow_v(this%phys,stat_info_sub)
        
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF ( num_colloid > 0 ) THEN
           
           CALL physics_get_colloid(this%phys,&
                colloids,stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Get current flow velocity.
        !----------------------------------------------------
        
        CALL statistic_compute_v_average(this%statis, &
             this%particles,stat_info_sub)
        CALL statistic_get_v_average(this%statis, &
             flow_v_current, stat_info_sub)
        
        !----------------------------------------------------
        ! Get current body force.
        !----------------------------------------------------
              
        CALL physics_get_body_force(this%phys,body_force, &
             stat_info_sub)
        
        !----------------------------------------------------
        ! Get in-/decrement of body force.
        !----------------------------------------------------
        
        Call physics_get_body_force_d(this%phys,&
             body_force_d,stat_info_sub)
        
        !--------------------------------------------------
        ! Dynamically changing the increament
        !(body_force_d) of body force.
        ! IF after first adjustment of body force,
        ! the flow velocity is smaller, but bigger after the
        ! second adjustment, body_force_d will be halfed.
        ! The similar rules applied to the case of 
        ! bigger than required at first adjustment,
        ! but smaller than required at second ajustment,
        ! body_force_d will be doubled.
        !--------------------------------------------------
        
        fd = physics_get_flow_direction(this%phys, &
             stat_info_sub)
              
        IF( flow_v_current(fd) < flow_v ) THEN
           
           body_force(1:num_dim) = &
                body_force(1:num_dim) + &
                body_force_d(1:num_dim)
           
        ELSE IF( flow_v_current(fd) > flow_v) THEN
           
           body_force(1:num_dim) = &
                body_force(1:num_dim) - &
                body_force_d(1:num_dim)
           
        END IF
        
        !----------------------------------------------------
        ! Set new body force.
        !----------------------------------------------------
        
        CALL physics_set_body_force(this%phys, &
             body_force,stat_info_sub)
        
        IF( num_colloid > 0 ) THEN
           
           CALL colloid_set_body_force(colloids, &
                body_force,stat_info_sub)
           
        END IF
        
        IF( stat_info_sub /=0 ) THEN
           PRINT *,"marching_adjust_flow_v : ", &
                "Adjusting body force failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
9999    CONTINUE
        
        
        IF(ASSOCIATED(body_force)) THEN
           DEALLOCATE(body_force)
        END IF
        
        IF(ASSOCIATED(body_force_d)) THEN
           DEALLOCATE(body_force_d)
        END IF
        
        IF(ASSOCIATED(flow_v_current)) THEN
           DEALLOCATE(flow_v_current)
        END IF
        
        RETURN
        
      END SUBROUTINE marching_adjust_flow_v
      
