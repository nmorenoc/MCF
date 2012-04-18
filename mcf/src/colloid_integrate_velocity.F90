      SUBROUTINE colloid_integrate_velocity(this,&
           step,dt,stat_info)
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
        ! Revision    : V0.3, 18.04.2012, 
        !               Adams-Bashforth method
        !               is implemented,
        !               inlcuding 1, 2, 3, 4, and 5 order.
        !
        !               V0.2 05.10.2009, including rotating
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
        INTEGER, INTENT(IN)             :: step
        REAL(MK), INTENT(IN)            :: dt        
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables:
        !----------------------------------------------------
        
        INTEGER                         :: itype, order, i

        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info  = 0
        itype = this%integrate_type

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

        
        IF ( this%translate ) THEN
        
           !-------------------------------------------------
           ! Save the previous translational velocity.
           !-------------------------------------------------

           i = itype
           
           DO WHILE ( i >= 2 )
           
              this%v(:,:,i)  =  this%v(:,:,i-1)
              
              i = i-1              
              
           END DO
           
           
           SELECT CASE (order)
              
           CASE (1)
              
              this%v(:,:,1) = &
                   this%v(:,:,1) + this%f(:,:,1) * dt
              
           CASE (2)
              
              this%v(:,:,1) = &
                   this%v(:,:,1) + &
                   ( 3.0_MK * this%f(:,:,1) - &
                   this%f(:,:,2)) * dt / 2.0_MK
              
           CASE (3)
              
              this%v(:,:,1) = &
                   this%v(:,:,1) + &
                   ( 23.0_MK * this%f(:,:,1) - &
                   16.0_MK * this%f(:,:,2) + &
                   5.0_MK * this%f(:,:,3) ) * dt / 12.0_MK
              
           CASE (4)
              
              this%v(:,:,1) = &
                   this%v(:,:,1) + &
                   ( 55.0_MK * this%f(:,:,1) - &
                   59.0_MK * this%f(:,:,2) + &
                   37.0_MK * this%f(:,:,3) - &
                   9.0_MK * this%f(:,:,4) ) * dt / 24.0_MK
              
           CASE (5) 
              
              this%v(:,:,1) = &
                   this%v(:,:,1) + &
                   ( 1901.0_MK * this%f(:,:,1) - &
                   2774.0_MK * this%f(:,:,2) + &
                   2616.0_MK * this%f(:,:,3) - &
                   1274.0_MK * this%f(:,:,4) + &
                   251.0_MK * this%f(:,:,5) ) * dt / 720.0_MK
              
           CASE DEFAULT
              
              PRINT *, "colloid_integrate_velocity: ",&
                   "integrator not available!"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! order
           
        END IF ! translate
        
        
        IF ( this%rotate ) THEN
           
           !-------------------------------------------------
           ! Save the previous rotational velocity.
           !-------------------------------------------------
           
           i = itype
           
           DO WHILE ( i >=  2 )
              
           this%omega(:,:,i)  =  this%omega(:,:,i-1)
           
           i = i -1
              
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
              
              PRINT *, "colloid_integrate_velocity: ",&
                   "integrator not available!"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! order
           
        END IF ! rotate
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_integrate_velocity
      
