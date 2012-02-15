      SUBROUTINE physics_init(this,d_ctrl,stat_info)
        !----------------------------------------------------
        ! Subroutine  :  physics_init
        !----------------------------------------------------
        !
        ! Purpose     : Construtor of Class Physics.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 15.07.2009, original version.
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
        
        TYPE(Physics),INTENT(OUT)               :: this
        TYPE(Control), INTENT(IN),TARGET        :: d_ctrl
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !----------------------------------------------------
        ! Init values are for flow around cylinder 
        ! with resolution 50*50
        ! dx(1:2) = 2.0e-3.
        ! This will be overritten, 
        ! if the physics.config file presents.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(this%ctrl)
        this%ctrl => d_ctrl
        
        this%num_species       = 2
        this%num_dim           = 2

        NULLIFY(this%min_phys)
        ALLOCATE(this%min_phys(2))
        this%min_phys(1:2) = 0.0_MK

        NULLIFY(this%max_phys)
        ALLOCATE(this%max_phys(2))
        this%max_phys(1)       = 1.0e-1_MK
        this%max_phys(2)       = 1.0e-1_MK
        
        NULLIFY(this%min_phys_t)
        ALLOCATE(this%min_phys_t(2))
        this%min_phys_t(1:2)   = this%min_phys(1:2)
        
        NULLIFY(this%max_phys_t)
        ALLOCATE(this%max_phys_t(2))
        this%max_phys_t(1:2)   = this%max_phys(1:2)
        
        this%lattice_type      = 1

        NULLIFY(this%num_part_dim)
        ALLOCATE(this%num_part_dim(2))
        
        this%num_part_dim(1)   = 50
        this%num_part_dim(2)   = 50
        
        NULLIFY(this%num_part_dim_t)
        ALLOCATE(this%num_part_dim_t(2))
        
        this%num_part_dim_t(1)   = 50
        this%num_part_dim_t(2)   = 50
        
        this%num_part_tot      = 2500
        
        NULLIFY(this%dx)
        ALLOCATE(this%dx(2))
        this%dx(1:2)           = &
             (this%max_phys(1:2)-this%min_phys(1:2)) / &
             this%num_part_dim(1:2)
        
        this%cut_off       = 9.6e-3_MK
        this%h             = 3.2e-3_MK
        this%dt            = -1.0_MK
        this%dt_c          = -1.0_MK
        this%dt_nu         = -1.0_MK
        this%fa_max        = -1.0_MK
        this%dt_f          = -1.0_MK
        this%step_start    = -1
        this%step_end      = -1
        this%step_current  = 0
        this%time_start    = 0.0_MK
        this%time_end      = 500.0_MK
        this%time_current  = 0.0_MK
        this%rho           = 1.0e+3_MK
        this%eta           = 1.0e-1_MK
        this%eta_coef      = 4.0_MK
        this%ksai          = 0.0_MK
        this%kt            = 0.0_MK
        this%c             = 2.0e-2_MK
        this%rho_ref       = this%rho*0.9_MK
        this%gamma         = 7.0_MK
        
        this%relax_type     = 0
        this%dt_relax       = 0.0_MK
        this%step_relax     = 0
        this%time_relax     = 0.0_MK
        this%disorder_level = 0.2_MK
        this%kt_relax       = 0.0_MK
        this%c_relax        = 0.0_MK

        this%tau  = 0.0_MK
	this%tau_sm = 0.0_MK
	this%k_sm = 0.0_MK
        this%n_p  = 0.0_MK
        this%kt_p = 0.0_MK
        this%eigen_dynamics    = .FALSE.
        NULLIFY(this%eval)
        this%eval_tolerance   = 0.0_MK        
        ALLOCATE(this%eval(4))
        NULLIFY(this%evec)
        ALLOCATE(this%evec(4,4))
        this%evec_normalize   = .FALSE.
        this%evec_tolerance   = 0.0_MK
        
        this%body_force_type = 1
        NULLIFY(this%body_force)
        ALLOCATE(this%body_force(2))
        this%body_force(1)    = 5.0e-5_MK
        this%body_force(2)    = 0.0_MK
        NULLIFY(this%body_force_d)
        ALLOCATE(this%body_force_d(2))      
        this%body_force_d(1:2)= 0.0_MK
        this%flow_direction   = 1
        this%flow_width       = 0.1_MK
        this%flow_v           = 0.0_MK
        this%flow_adjust_freq = 0
      
        this%num_colloid   = 0
        NULLIFY(this%colloids)
        
        NULLIFY(this%bcdef)
        ALLOCATE(this%bcdef(4))
        this%bcdef(1:4)       = ppm_param_bcdef_periodic
        
        NULLIFY(this%boundary)

        
        RETURN
          
      END SUBROUTINE physics_init
     
      
      SUBROUTINE physics_display_parameters(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : physics_display parameters
        !----------------------------------------------------
        !
        ! Purpose     : To display physics parameters.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 15.07.2009, original version.
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

        TYPE(Physics),INTENT(IN)      :: this
        INTEGER,INTENT(OUT)           :: stat_info
        
        
        LOGICAL                       :: relax_run
        LOGICAL                       :: flow_v_fixed
        LOGICAL                       :: Newtonian
        
        INTEGER                       :: j,num_dim
        INTEGER                       :: stat_info_sub
          
        stat_info     = 0
        stat_info_sub = 0
        
        relax_run    = &
             control_get_relax_run(this%ctrl,stat_info_sub)
        flow_v_fixed = &
             control_get_flow_v_fixed(this%ctrl,stat_info_sub)      
        Newtonian = &
             control_get_Newtonian(this%ctrl,stat_info_sub)
        num_dim   = this%num_dim
        
        PRINT *, '------------------Start------------------'
        PRINT *, '     Physics parameters'
        PRINT *, '-----------------------------------------'
        
        PRINT *, "num_species      : ", this%num_species
        PRINT *, "num_dim          : ", num_dim
        PRINT *, "min_phys         : "
        PRINT *, this%min_phys(1:num_dim)
        PRINT *, "max_phys         : "
        PRINT *, this%max_phys(1:num_dim)
        PRINT *, "min_phys_t       : "
        PRINT *, this%min_phys_t(1:num_dim)
        PRINT *, "max_phys_t       : "
        PRINT *, this%max_phys_t(1:num_dim)
      
        IF(this%num_dim ==2) THEN
           
           SELECT CASE (this%lattice_type)
           CASE (1)
              PRINT *, "lattice type     : ", "  Square Lattice."
           CASE (2)
              PRINT *, "lattice type     : ", "  Staggered Lattice."
           CASE (3)
              PRINT *, "lattice type     : ", "  Hexagonal Lattice."
           END SELECT
           
        ELSE IF(this%num_dim == 3) THEN
           
           SELECT CASE (this%lattice_type)
           CASE (1)
              PRINT *, "lattice type     : ", "  Simple Cubic Lattice."
           CASE (2)
              PRINT *, "lattice type     : ", "  Body Center Lattice."
           CASE (3)
              PRINT *, "lattice type     : ", "  Face Center Lattice."       
           END SELECT
           
        END IF
        
        PRINT *, "num_part_dim     : "
        PRINT *, this%num_part_dim(1:num_dim)
        PRINT *, "num_part_dim_t   : "
        PRINT *, this%num_part_dim_t(1:num_dim)
        PRINT *, "num_part_tot     : ", this%num_part_tot
        PRINT *, "dx               : "
        PRINT *, this%dx(1:num_dim)
        PRINT *, "cut_off          : ", this%cut_off
        PRINT *, "smoothing length : ", this%h
        PRINT *, "dt               : ", this%dt
        PRINT *, "dt_c             : ", this%dt_c
        PRINT *, "dt_nu            : ", this%dt_nu
        PRINT *, "dt_f             : ", this%dt_f
        PRINT *, "step_start       : ", this%step_start
        PRINT *, "step_end         : ", this%step_end     
        PRINT *, "time_start       : ", this%time_start
        PRINT *, "time_end         : ", this%time_end
        PRINT *, "rho              : ", this%rho
        PRINT *, "eta              : ", this%eta
        PRINT *, "eta_coef         : ", this%eta_coef
        PRINT *, "ksai             : ", this%ksai
        PRINT *, "kt               : ", this%kt
        PRINT *, "c                : ", this%c
        PRINT *, "rho_ref          : ", this%rho_ref
        PRINT *, "gamma            : ", this%gamma
	PRINT *, "tau_sm           : ", this%tau_sm
	PRINT *, "k_sm             : ", this%k_sm
        
        IF ( relax_run ) THEN

           PRINT *, "relax_type       : ", this%relax_type
           PRINT *, "dt_relax         : ", this%dt_relax
           PRINT *, "dt_c_relax       : ", this%dt_c_relax
           PRINT *, "step_relax       : ", this%step_relax
           PRINT *, "time_relax       : ", this%time_relax
           PRINT *, "disorder_level   : ", this%disorder_level
           PRINT *, "kt_realx         : ", this%kt_relax
           PRINT *, "c_realx          : ", this%c_relax

        END IF
        
        IF (.NOT. Newtonian ) THEN
           
           PRINT *, "tau              : ", this%tau
           PRINT *, "n_p              : ", this%n_p
           PRINT *, "kt_p             : ", this%kt_p
           PRINT *, "eigen_dynamics   : ", this%eigen_dynamics
           IF( this%eigen_dynamics) THEN
              PRINT *, "eigenvalues      : "
              PRINT *, this%eval(1:num_dim)
              PRINT *, "eval_tolerance   : ", this%eval_tolerance
              PRINT *, "eigenvectors     : "
              DO j = 1, num_dim
                 PRINT *, this%evec(1:num_dim,j)
              END DO
              PRINT *, "evec_normalize   : ", this%evec_normalize
              IF ( this%evec_normalize ) THEN
                 PRINT *, "evec_tolerance : ", this%evec_tolerance
              END IF
           END IF
           
        END IF
        
        PRINT *, "body_force_type  : ", this%body_force_type
        PRINT *, "body_force       : "
        PRINT *, this%body_force(1:num_dim)
        
        IF ( flow_v_fixed ) THEN
           
           PRINT *, "body_force_d     : "
           PRINT *,  this%body_force_d(1:num_dim)
           
           PRINT *, "flow_direction   : ", this%flow_direction
           PRINT *, "flow_width       : ", this%flow_width
           PRINT *, "flow_v           : ", this%flow_v
           PRINT *, "flow_adjust_freq : ", this%flow_adjust_freq
           
        END IF
        
        PRINT *, "num_colloid      : ", this%num_colloid
        
        IF (this%num_colloid > 0) THEN
           
           CALL colloid_display_parameters(this%colloids,stat_info_sub)
           
        END IF
        
        CALL boundary_display_parameters(this%boundary,stat_info_sub)

        PRINT *, '-------------------End-------------------'
        
        RETURN
        
      END SUBROUTINE physics_display_parameters
      
      
