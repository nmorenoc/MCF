      SUBROUTINE colloid_compute_rotation_vector(this,&
           dt,accuracy_order,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_rotation_vector
        !----------------------------------------------------
        !
        ! Purpose     : Compute the rotation vector
        !
        ! Remark      : Colloid are modeled as rigid body.
        !
        ! Reference  : Chen et. al. 2006, physics of fluids
        !              wikipedia
        !
        ! Revision   : V0.1  14.10.2009, original.
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
        INTEGER                         :: dim,i
        REAL(MK),DIMENSION(3)           :: axis
        REAL(MK)                        :: phi
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        dim           = this%num_dim

#if __POSITION_FIXED
#else        
        IF ( this%rotate ) THEN

           !-------------------------------------------------
           ! Select different accuracy oder:
           ! 1 velocity contribution;
           ! 2 velocity + acceleration contribution.
           !-------------------------------------------------
           
           SELECT CASE (accuracy_order)
              
           CASE(1)
           
              DO i = 1, this%num_colloid
                 
                 !-------------------------------------------
                 ! Calculate the roation vector at this time
                 ! step and normalize it.
                 !-------------------------------------------

                 axis(1:3) = this%omega(1:3,i) * dt
                 
                 phi = SQRT(DOT_PRODUCT(axis(1:3),axis(1:3)))
                 
                 IF ( phi < mcf_machine_zero ) THEN
                    phi       = 0.0_MK
                    axis(1)   = 1.0_MK
                    axis(2:3) = 0.0_MK
                 ELSE
                    axis(1:3) = axis(1:3) / phi
                 END IF
                 
                 this%rot_vector(1:3,i) = axis(1:3)
                 this%rot_vector(4,i)   = phi
                 
                 !-------------------------------------------
                 ! Only usefull for 2D,
                 ! since 3D roation angle can not be simply
                 ! added up for accumulating rotated angle.
                 !-------------------------------------------
                 
                 this%theta(3,i) = this%theta(3,i) + phi
                 
              END DO
              
           CASE(2)
              
              DO i = 1, this%num_colloid
                 
                 !-------------------------------------------
                 ! Calculate the roation vector at this time
                 ! step.
                 !-------------------------------------------
                 
                 axis(1:3) = this%omega(1:3,i) * dt + &
                      0.5_MK * this%alpha(1:3,i) * dt**2

                 phi = SQRT(DOT_PRODUCT(axis(1:3),axis(1:3)))
                 
                 IF ( phi < mcf_machine_zero ) THEN
                    phi   = 0.0_MK
                    axis(1)   = 1.0_MK
                    axis(2:3) = 0.0_MK
                 ELSE
                    axis(1:3) = axis(1:3) / phi
                 END IF
                 
                 this%rot_vector(1:3,i) = axis(1:3)
                 this%rot_vector(4,i)   = phi
                 
                 !-------------------------------------------
                 ! Only usefull for 2D, 
                 ! since 3D roation vector can not be simply
                 ! added up for accumulating rotated angle.
                 !-------------------------------------------
                 
                 this%theta(3,i) = this%theta(3,i) + phi
                 
              END DO ! i = 1, num_colloid
              
           END SELECT ! accuracy_order
           
        END IF ! rotate
        
#endif
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_compute_rotation_vector
      
