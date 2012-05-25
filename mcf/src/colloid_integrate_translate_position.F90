      SUBROUTINE colloid_integrate_translate_position(this,&
           step,dt,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_integrate_translate_position
        !----------------------------------------------------
        !
        ! Purpose     : Integrate the positions of centers of
        !               colloids.
        !
        ! Remark      : Colloid are modeled as rigid body.
        !
        !
        !  Revision   : V0.2, 18.04.2012, 
        !               Adams-Bashforth method is 
        !               implemented inlcuding 
        !               1, 2, 3, 4, and 5 order accuracy.
        !
        !               V0.1  23.06.2009, original.
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
        
        INTEGER                         :: itype, order
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        itype     = this%integrate_AB
        
#ifdef __PARTICLES_POSITION_FIXED
#else
        
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
        
        IF( this%translate ) THEN
           
           SELECT CASE (order)
              
           CASE (1)
              
              this%x(:,:) = &
                   this%x(:,:) + this%v(:,:,1) * dt
              
           CASE (2)
              
              this%x(:,:) = &
                   this%x(:,:) + &
                   ( 3.0_MK * this%v(:,:,1) - &
                   this%v(:,:,2) ) * dt / 2.0_MK
              
           CASE (3)
              
              this%x(:,:) = &
                   this%x(:,:) + &
                   ( 23.0_MK * this%v(:,:,1) - &
                   16.0_MK * this%v(:,:,2) + &
                   5.0_MK * this%v(:,:,3) ) * dt / 12.0_MK
              
           CASE (4)
              
              this%x(:,:) = &
                   this%x(:,:) + &
                   ( 55.0_MK * this%v(:,:,1) - &
                   59.0_MK * this%v(:,:,2) + &
                   37.0_MK * this%v(:,:,3) - &
                   9.0_MK * this%v(:,:,4) ) * dt / 24.0_MK
              
           CASE (5) 
              
              this%x(:,:) = &
                   this%x(:,:) + &
                   ( 1901.0_MK * this%v(:,:,1) - &
                   2774.0_MK * this%v(:,:,2) + &
                   2616.0_MK * this%v(:,:,3) - &
                   1274.0_MK * this%v(:,:,4) + &
                   251.0_MK * this%v(:,:,5) ) * dt / 720.0_MK
              
           CASE DEFAULT
              
              PRINT *, "colloid_integrate_position: ",&
                   "integrator not available!"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! order
           
        END IF ! translate
#endif
   
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_integrate_translate_position
      
