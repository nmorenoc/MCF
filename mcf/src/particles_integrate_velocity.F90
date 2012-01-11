      SUBROUTINE particles_integrate_velocity(this,&
           num,dt,lambda,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_integrate_velocity
        !----------------------------------------------------
        !
        ! Purpose     : Integrate velocity of particles
        !               using accleration f with lamda, 
        !               which is a prediction or
        !               correction coefficient.
        !                  
        !                  
        ! Reference   :
        !
        ! Remark      : colloidal boundary particle may rotate,
        !               needs to be done seperately.
        !               wall boundary particles may have
        !               different behavior, such as,
        !               oscillating shear, needs to be done
        !               seperately.
        !
        ! Revisions   : V0.4 22.11.2010, integrate velocity of
        !               only fluid particles.
        ! 
        !               V0.3 24.08.2010, integrate velocity of
        !               particles, except colloidal boundary
        !               particles.
        !
        !               V0.2 09.07.2009, 
        !               check again the work flow is correct and
        !               supply with more comments for code.
        !
        !               V0.1 01.04.2009, original version.
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
        !  Arguments
        !
        !  this       : an object of Particles Class.
        !  num        : first num of particles needed to be updated.
        !  dt         : time step.
        !  lamda      : coefficient of acceleration.
        !  stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        REAL(MK), INTENT(IN)                    :: dt
        REAL(MK), INTENT(IN)                    :: lambda
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        !  Local variables
	!----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim,i
   

        !----------------------------------------------------
        !  Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        
        dim  = physics_get_num_dim(this%phys,stat_info_sub)
        
        !----------------------------------------------------
        ! Update velcoity of fluid particles.
        !----------------------------------------------------
        
        DO i = 1, num
           
           IF ( this%id(this%sid_idx,i) == mcf_particle_type_fluid ) THEN
              
              this%v(1:dim,i) = &
                   this%v(1:dim,i) + &
                   lambda * this%f(1:dim,i) * dt
              
           END IF
           
        END DO ! i =1, num

9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE particles_integrate_velocity
      
