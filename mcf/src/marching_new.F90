      SUBROUTINE marching_init(this,d_io,d_particles,stat_info)
        !----------------------------------------------------
        ! Subroutine  : marching_init
        !----------------------------------------------------
        !
        ! Purpose     : Constructor of Marching.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 15.07 2009, original version.
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
        ! Arguments :
        !
        ! this       : an object of Marching Class.
        ! stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Marching), INTENT(OUT)             :: this
        TYPE(IO), INTENT(IN), TARGET            :: d_io
        TYPE(Particles), INTENT(IN), TARGET     :: d_particles
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        LOGICAL                                 :: flow_v_fixed
        LOGICAL                                 :: p_energy
        INTEGER                                 :: num_dim
        
        !----------------------------------------------------
        ! Initializatoin of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        
        CALL particles_get_ctrl(d_particles,this%ctrl,stat_info_sub)
        CALL particles_get_phys(d_particles,this%phys,stat_info_sub)
        this%io             => d_io
        this%particles      => d_particles
        CALL particles_get_tech(this%particles,this%tech,stat_info_sub)
        this%integrate_type = &
             control_get_integrate_type(this%ctrl,stat_info_sub)

        
        flow_v_fixed = &
             control_get_flow_v_fixed(this%ctrl,stat_info_sub)
        p_energy = &
             control_get_p_energy(this%ctrl,stat_info_sub)
        num_dim = physics_get_num_dim(this%phys,stat_info_sub)
        
        CALL statistic_new(this%statis, &
             this%ctrl,num_dim, stat_info_sub)
        
        RETURN
        
      END SUBROUTINE marching_init


      
      SUBROUTINE marching_display_parameters(this,stat_info)
        
        TYPE(Marching), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        

        
        stat_info     = 0
        
        
        PRINT *, '------------------Start------------------'
        PRINT *,'      Marching Parameters '
        PRINT *, '-----------------------------------------'
        
        SELECT CASE (this%integrate_type)
           
        CASE (1)
           
           PRINT *, "integration type : ", &
                "Explicit Euler"
           
        CASE (2)
           
           PRINT *, "integration type : ", &
                "Modified velocity Verlet"
                     
        END SELECT
        
        PRINT *, '-------------------End-------------------'
        
        RETURN          
        
      END SUBROUTINE marching_display_parameters
      
      
