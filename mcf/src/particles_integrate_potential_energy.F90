      SUBROUTINE particles_integrate_potential_energy(this,&
           num,dt,lamda,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_integrate_potential_energy
        !----------------------------------------------------
        !
        ! Purpose     : Integrate potential energy of first 
        !               num particles, using accleration au
        !               with a coefficient lamda.
        !                  
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.2 09.07.2009, origianl version.
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
        ! Arguments
        !
        ! this        : an object of Particles Class.
        ! num         : first num of particles needed to be updated.
        ! lamda       : coefficient of accleartion.
        ! stat_info   : return flag of status.
        !----------------------------------------------------
       
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        REAL(MK), INTENT(IN)                    :: dt
        REAL(MK), INTENT(IN)                    :: lamda
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        

        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        
        dim  = physics_get_num_dim(this%phys,stat_info_sub)
        
        
        !----------------------------------------------------
        ! Integrate particles' potential energy.
        !----------------------------------------------------
        
        this%u(1:num) = this%u(1:num) + &
             this%au(1:num) * dt * lamda
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE particles_integrate_potential_energy
      
