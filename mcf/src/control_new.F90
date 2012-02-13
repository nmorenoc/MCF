      SUBROUTINE control_init_default(this,stat_info)
        !----------------------------------------------------
        !  Subroutine   :  control_init
        !----------------------------------------------------
        !
        !  Purpose      : Construtor of Class Control.
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : V0.1 15.07.2009, original version.
        !
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        TYPE(Control),INTENT(OUT)       :: this
        INTEGER,INTENT(OUT)             :: stat_info
        
        !----------------------------------------------------
        ! Init values are for flow around cylinder
        ! using
        ! quintic spline
        ! non-symmetry force calculation
        ! Morris's formulation of density and force
        ! Euler integration.
        !----------------------------------------------------
        
        stat_info = 0
        
        this%job_name            = ""
        this%job_submit_date     = ""
        CALL DATE_AND_TIME(this%job_submit_date)
        this%job_execute_date     = ""
        this%job_execute_time     = ""
        this%job_execute_zone     = ""
        CALL DATE_AND_TIME(this%job_execute_date, &
             this%job_execute_time, this%job_execute_zone)
        this%debug_flag          = 1
        this%relax_run           = .FALSE.
        this%read_external       = .FALSE.
        this%kernel_type         = 1
        this%symmetry            = .FALSE.
        this%rhs_density_type    = 1
        this%dynamic_density_ref = .FALSE.
        this%stateEquation_type  = 1
        this%Newtonian           = .TRUE.
        this%Brownian            = .FALSE.
        this%random_seed         = -1        
        this%rhs_force_type      = 2
        this%pp_interact_cc      = .FALSE.
        this%pp_interact_cw      = .FALSE.
        this%cc_lub_type         = 0
        this%cc_repul_type       = 0
        this%cw_lub_type         = 0
        this%cw_repul_type       = 0
        this%stress_tensor       = .FALSE.
        this%p_energy            = .FALSE.
        this%flow_v_fixed        = .FALSE.
        this%integrate_type      = 1
        this%adaptive_dt         = 0
        this%write_output        = 1
        this%write_restart       = 0

        RETURN
          
      END SUBROUTINE control_init_default
     
      
      SUBROUTINE control_display_parameters(this,stat_info)
        !----------------------------------------------------
        !  Subroutine   :  control_init
        !----------------------------------------------------
        !
        !  Purpose      : To display control parameters
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : V0.1 15.07.2009, original version.
        !
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        TYPE(Control),INTENT(IN)        :: this
        INTEGER,INTENT(OUT)             :: stat_info
        
        stat_info = 0
        
        PRINT *, '------------------Start------------------'
        PRINT *, '     Control  parameters'
        PRINT *, '-----------------------------------------'
        
        PRINT *, "job_name           : ", &
             this%job_name(1:LEN_TRIM(this%job_name))
        PRINT *, "job_submit_date    : ", &
             this%job_submit_date(1:LEN_TRIM(this%job_submit_date))
        PRINT *, "job_execute_date   : ", &
             this%job_execute_date(1:LEN_TRIM(this%job_execute_date))
        PRINT *, "job_execute_time   : ", &
             this%job_execute_time(1:LEN_TRIM(this%job_execute_time))
        PRINT *, "job_execute_zone   : ", &
             this%job_execute_zone(1:LEN_TRIM(this%job_execute_zone))
     
        PRINT *, "debug_flag         : ", this%debug_flag
        
        PRINT *, "relax run          : ", this%relax_run
        
        PRINT *, "read external      : ", this%read_external
        
        SELECT CASE(this%kernel_type)
        CASE (1)
           PRINT *, "kernel_type        : ", "Quintic Spline"
        CASE (2)
           PRINT *, "kernel_type        : ", "Lucy kernel"
        CASE DEFAULT
           PRINT *, "kernel_type        : ", &
                this%kernel_type, " not available"
           stat_info = -1
           GOTO 9999
        END SELECT
        
        PRINT *, "Inter-symmetry     : ", this%symmetry
        
        SELECT CASE(this%rhs_density_type)
        CASE (1)
           PRINT *, "rhs_density_type   : ", &
                "Summation of mass density"
        CASE (2)
           PRINT *, "rhs_density_type   : ", &
                "Summation of number density"
        CASE DEFAULT
           PRINT *, "rhs_density_type   : ", &
                this%rhs_density_type, " not available"
           stat_info = -1
           GOTO 9999
        END SELECT
        
        PRINT *, "dynamic_density_ref: ", this%dynamic_density_ref
        
        SELECT CASE(this%stateEquation_type)
        CASE (1)
           PRINT *, "stateEquation_type : ", "Morris et al., 1997"
        CASE (2)
           PRINT *, "stateEquation_type : ", "Batchelor. 1967"
        CASE DEFAULT
           PRINT *, "stateEquation_type : ", &
                this%stateEquation_type, " not available"
        END SELECT
  
        PRINT *, "Newtonian fluid    : ", this%Newtonian
        
        PRINT *, "Brownian motion    : ", this%Brownian
        PRINT *, "random seed        : ", this%random_seed

        SELECT CASE(this%rhs_force_type)
        CASE (1)
           PRINT *, "rhs_force_type     : ", &
                "Morris et al., J. Comput. Phys. 1997"
        CASE (2)
           PRINT *, "rhs_force_type     : ", &
                "Espanol and Revenga, Phys. Rev. E 2003"
        CASE (3)
           PRINT *, "rhs_force_type     : ", &
                "Hu and Adams, J. Comput. Phys. 2006"
        CASE (4)
           PRINT *, "rhs_force_type     : ", &
                "Hu and Adams, Phys. Fluids 2006"
        CASE DEFAULT
           PRINT *, "rhs_force_type     : ", &
                this%rhs_force_type, " not available"           
        END SELECT

        PRINT *, "pp_interact_cc     : ", this%pp_interact_cc
        
        PRINT *, "pp_interact_cw     : ", this%pp_interact_cw
        
        PRINT *, "cc_lub_type        : ", this%cc_lub_type
        
        PRINT *, "cc_repul_type      : ", this%cc_repul_type
        
        PRINT *, "cw_lub_type        : ", this%cw_lub_type
        
        PRINT *, "cw_repul_type      : ", this%cw_repul_type
        
        PRINT *, "stress tensor      : ", this%stress_tensor
        
        PRINT *, "potential energy   : ", this%p_energy
        
        PRINT *, "flow velocity fixed: ", this%flow_v_fixed
        
        SELECT CASE(this%integrate_type)
        CASE (1)
           PRINT *, "integrate_type     : ", &
                "Explicit Euler"
        CASE (2)
           PRINT *, "integrate_type     : ", &
                "Modified velocity Verlet"
        CASE DEFAULT
           PRINT *, "integrate_type     : ", &
                this%integrate_type, " not available"
        END SELECT
        
        PRINT *, "adaptive dt        : ", this%adaptive_dt
        
        PRINT *, "write output       : ", this%write_output
        
        PRINT *, "write restart      : ", this%write_restart
        
        
        PRINT *, '-------------------End-------------------'
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE control_display_parameters
