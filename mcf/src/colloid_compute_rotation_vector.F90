      SUBROUTINE colloid_compute_rotation_vector(this,&
           step,dt,stat_info)
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
        INTEGER, INTENT(IN)             :: step
        REAL(MK), INTENT(IN)            :: dt
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim,i
        INTEGER                         :: itype, order
        REAL(MK),DIMENSION(3)           :: axis
        REAL(MK)                        :: phi
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        dim           = this%num_dim
        itype         = this%integrate_type
        
#if __PARTICLES_POSITION_FIXED
#else        

        !----------------------------------------------------
        ! Select different accuracy oder for
        ! When the step is smaller then integrator order,
        ! a lower order (step) integrator is used.
        ! When the step is bigger than or equal to integrator
        ! order, the actual order (itype) is used.
        !----------------------------------------------------
    
        IF ( step < itype ) THEN
           order = step
        ELSE
           order = itype             
        END IF

        IF ( this%rotate ) THEN
           
           DO i = 1, this%num_colloid
              
              !----------------------------------------------
              ! Calculate the roation axis at this time
              !----------------------------------------------
              
              SELECT CASE (order)
                 
              CASE (1)
                 
                 axis(1:3) = this%omega(1:3,i,1) * dt
                 
              CASE (2)   
                 
                 axis(1:3) = ( 3.0_MK * this%omega(1:3,i,1) - &
                      this%omega(1:3,i,2)) * dt / 2.0_MK
                 
              CASE (3)
                 
                 axis(1:3) = ( 23.0_MK * this%omega(1:3,i,1) - &
                      16.0_MK * this%omega(1:3,i,2) + &
                      5.0_MK * this%omega(1:3,i,3) ) * dt / 12.0_MK
                 
              CASE (4)
                 
                 axis(1:3) =  ( 55.0_MK * this%omega(1:3,i,1) - &
                      59.0_MK * this%omega(1:3,i,2) + &
                      37.0_MK * this%omega(1:3,i,3) - &
                      9.0_MK * this%omega(1:3,i,4) ) * dt / 24.0_MK
                 
              CASE (5)
                 
                 axis(1:3) = ( 1901.0_MK* this%omega(1:3,i,1) - &
                      2774.0_MK * this%omega(1:3,i,2) + &
                      2616.0_MK * this%omega(1:3,i,3) - &
                      1274.0_MK * this%omega(1:3,i,4) + &
                      251.0_MK * this%omega(1:3,i,5) ) * dt / 720.0_MK
                 
              CASE DEFAULT
                 
                 PRINT *, "colloid_compute_rotation_vector: ",&
                      "integrator not available!"
                 stat_info = -1
                 GOTO 9999
                 
              END SELECT ! order
              
              !-----------------------------------------------
              ! Calculate the roation vector at this time
              ! step and normalize it.
              !----------------------------------------------
              
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
              
           END DO ! num_colloid
           
        END IF ! rotate
        
#endif
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_compute_rotation_vector
      
#if 0      
      ! to delete
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
#endif
