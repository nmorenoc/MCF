      SUBROUTINE colloid_integrate_rotate_velocity(this,&
           step,dt,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_integrate_rotate_velocity
        !----------------------------------------------------
        !
        ! Purpose     : Integrate rotating velocity of 
        !               colloids.
        !
        ! Remark      : Colloid are modelled as rigid body.
        !               
        !
        ! Revision    : V0.3, 18.04.2012, 
        !               Adams-Bashforth method is 
        !               implemented inlcuding 
        !               1, 2, 3, 4, and 5 order accuracy.
        !
        !               V0.2 05.10.2009, including rotating
        !               velocity.(deleted 1.1.2012)
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
        INTEGER, INTENT(IN)             :: step
        REAL(MK), INTENT(IN)            :: dt        
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables:
        !----------------------------------------------------
        
        INTEGER                         :: itype, order, i

        !----------------------------------------------------
        ! Initialization of variables:
        !----------------------------------------------------
        
        stat_info = 0
        itype     = this%integrate_AB

        !----------------------------------------------------
        ! Select different accuracy oder:
        ! when the step is smaller then desired accuracy order,
        ! a lower order (step) integrator is used, as we have
        ! no more information about history.
        ! For example, at 1st step using explicit Euler,
        !
        ! When the step is bigger than or equal to desired
        ! accuracy order, the actual deisired order (itype) 
        ! can be used and will be used.
        !----------------------------------------------------
           
        IF ( step < itype ) THEN
           order = step
        ELSE
           order = itype             
        END IF

        IF ( this%rotate ) THEN
           
           !-------------------------------------------------
           ! Save the rotational velocity at previous 
           ! time steps.
           !-------------------------------------------------
           
           i = itype
           
           DO WHILE ( i >=  2 )
              
              this%omega(:,:,i)  =  this%omega(:,:,i-1)
           
              i = i - 1
              
           END DO
           
           SELECT CASE (order)
              
           CASE (1)
              
              this%omega(:,:,1) = &
                   this%omega(:,:,1) + this%alpha(:,:,1) * dt
              
           CASE (2)
              
              this%omega(:,:,1) = &
                   this%omega(:,:,1) + &
                   ( 3.0_MK * this%alpha(:,:,1) - &
                   this%alpha(:,:,2) ) * dt / 2.0_MK
              
           CASE (3)
              
              this%omega(:,:,1) = &
                   this%omega(:,:,1) + &
                   ( 23.0_MK * this%alpha(:,:,1) - &
                   16.0_MK * this%alpha(:,:,2) + &
                   5.0_MK * this%alpha(:,:,3) ) * dt / 12.0_MK
              
           CASE (4)
              
              this%omega(:,:,1) = &
                   this%omega(:,:,1) + &
                   ( 55.0_MK * this%alpha(:,:,1) - &
                   59.0_MK * this%alpha(:,:,2) + &
                   37.0_MK * this%alpha(:,:,3) - &
                   9.0_MK * this%alpha(:,:,4) ) * dt / 24.0_MK
              
           CASE (5) 
              
              this%omega(:,:,1) = &
                   this%omega(:,:,1) + &
                   ( 1901.0_MK * this%alpha(:,:,1) - &
                   2774.0_MK * this%alpha(:,:,2) + &
                   2616.0_MK * this%alpha(:,:,3) - &
                   1274.0_MK * this%alpha(:,:,4) + &
                   251.0_MK * this%alpha(:,:,5) ) * dt / 720.0_MK
              
           CASE DEFAULT
              
              PRINT *, "colloid_integrate_rotate_velocity: ",&
                   "integrator not available!"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! order
           
        END IF ! rotate
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_integrate_rotate_velocity
      
