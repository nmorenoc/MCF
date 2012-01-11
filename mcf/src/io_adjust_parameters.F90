      SUBROUTINE io_adjust_parameters(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_adjust_parameters
        !----------------------------------------------------
        !
        ! Purpose     : Adjust io parameters after
        !               they have values, and resolve
        !               the small conflicts.
        !      
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 23.07.2009, original version.
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
        !----------------------------------------------------
        
        TYPE(IO), INTENT(INOUT)         :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local parameters
        !----------------------------------------------------
        INTEGER                         :: stat_info_sub

        
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        
        SELECT CASE (this%write_output) 
           
        CASE (1)
           this%output_particles_freq_time    = -1.0_MK
           this%output_conformation_freq_time = -1.0_MK
           this%statistic_freq_time = -1.0_MK
           this%boundary_freq_time  = -1.0_MK
           this%colloid_freq_time   = -1.0_MK
        CASE (2)  
           this%output_particles_freq_step    = -1 
           this%output_conformation_freq_step = -1
           this%statistic_freq_step  = -1
           this%boundary_freq_step   = -1
           this%colloid_freq_step    = -1                   
        END SELECT
        
        this%step_start = &
             physics_get_step_start(this%phys,stat_info_sub)
        this%time_start = &
             physics_get_time_start(this%phys,stat_info_sub)
      
        SELECT CASE ( this%write_restart)
        CASE (1)
           this%restart_freq_time      = -1.0_MK
           this%restart_freq_time_wall = -1.0_MK
        CASE (2)
           this%restart_freq_step       = -1
           this%restart_freq_time_wall  = -1.0_MK
        CASE (3)
           this%restart_freq_step  = -1
           this%restart_freq_time  = -1.0_MK
        END SELECT
           
        RETURN
        
      END SUBROUTINE io_adjust_parameters
      
