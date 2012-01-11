      LOGICAL FUNCTION physics_check_parameters(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : physics_check_parameters
        !----------------------------------------------------
        !
        ! Purpose     : Check if physics parameters are
        !               resonable.
        !
        ! Reference   :
        !
        ! Remark      : need to be extended to 
        !               check all parameters.
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
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: lcheck
        LOGICAL                         :: relax_run
        LOGICAL                         :: read_external
        LOGICAL                         :: symmetry
        LOGICAL                         :: Newtonian
        LOGICAL                         :: Brownian
        INTEGER                         :: adaptive_dt
        INTEGER                         :: dim
        INTEGER                         :: num_peri
        INTEGER                         :: num_sym
        INTEGER                         :: num_wall_sym
        INTEGER                         :: num_wall_solid
        INTEGER                         :: num_wall
        INTEGER                         :: num_le
        TYPE(Colloid), POINTER          :: colloids
        INTEGER                         :: i
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        lcheck    = .TRUE.        
        
        relax_run     = &
             control_get_relax_run(this%ctrl,stat_info_sub)
        read_external = &
             control_get_read_external(this%ctrl,stat_info_sub)
        symmetry      = &
             control_get_symmetry(this%ctrl,stat_info_sub)
        Newtonian      = &
             control_get_Newtonian(this%ctrl,stat_info_sub)
        Brownian      = &
             control_get_Brownian(this%ctrl,stat_info_sub)
        adaptive_dt   = &
             control_get_adaptive_dt(this%ctrl,stat_info_sub)
     
        dim  = this%num_dim
        num_peri  = &
             boundary_get_num_peri(this%boundary,stat_info_sub)
        num_sym   = &
             boundary_get_num_sym(this%boundary,stat_info_sub)
        num_wall_sym = &
             boundary_get_num_wall_sym(this%boundary,stat_info_sub)
        num_wall_solid = &
             boundary_get_num_wall_solid(this%boundary,stat_info_sub)
        num_wall = &
             boundary_get_num_wall(this%boundary,stat_info_sub)
        num_le   = &
             boundary_get_num_le(this%boundary,stat_info_sub)
       
        NULLIFY(colloids)
        
        !----------------------------------------------------
        ! Check if the number of species is resonable.
        !----------------------------------------------------
        
        IF ( this%num_species < 1 .OR. &
             this%num_species > 2) THEN
           PRINT *, "phsyics_check_parameters : ", &
                "number of species can only be 1 or 2 !"
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check if the number of dimension is resonable.
        !----------------------------------------------------
        
        IF ( dim < 2 .OR. dim > 3) THEN
           PRINT *, "phsyics_check_parameters : ", &
                "number of dimension can only be 2 or 3 !"
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! start step can not be bigger than end step.
        !----------------------------------------------------
        
        IF ( this%step_start > this%step_end ) THEN
           PRINT *, "phsyics_check_parameters : ", &
                "step_step > step_end, wrong !"
           lcheck = .FALSE.
           GOTO 9999
        END IF

        !----------------------------------------------------
        ! start time can not be bigger than end time.
        !----------------------------------------------------
        
        IF ( this%time_start > this%time_end ) THEN
           PRINT *, "phsyics_check_parameters : ", &
                "time_step > time_end, wrong !"
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! At least one of step_* couple or 
        ! time_* couple should be given as
        ! starting and ending point for simulation.
        !----------------------------------------------------
        
        IF ( this%step_end < 0 .AND. &
             this%time_end < 0.0_MK ) THEN
           PRINT *, "phsyics_check_parameters : ", &
                "One of step_end and time_end must be non-negative !"
           lcheck = .FALSE.
           GOTO 9999
        END IF
        

        IF ( (read_external .AND. adaptive_dt>0 ) &
             .OR. (.NOT. read_external) ) THEN
           
           !----------------------------------------------------
           ! For reading external and using adaptive dt,
           ! or for non-reading external, 
           ! check step and time,
           ! because only one pair, i.e, either step paris
           ! or time pair, has be given.
           ! For reading external run, this doesn't matter.
           !----------------------------------------------------

           IF ( this%step_end >= 0 .AND. &
                this%time_end >= 0.0_MK ) THEN
              PRINT *, "phsyics_check_parameters : ", &
                   "Only one of step and time should be non-negative !"
              lcheck = .FALSE.
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! IF step is given from input,
           ! then 0<=step_start<= step_end.
           !-------------------------------------------------
           
           IF( this%step_end >= 0 ) THEN
              
              IF ( this%step_start < 0 ) THEN
                 
                 PRINT *, "phsyics_check_parameters : ", &
                      " step_start does not match step_end !"
                 lcheck = .FALSE.
                 GOTO 9999
                 
              END IF

           END IF ! step_end >= 0
           
           !-------------------------------------------------
           ! IF time is given from input,
           ! then 0<=time_start<= time_end.
           !-------------------------------------------------
           
           IF ( this%time_end >=0.0_MK ) THEN
              
              IF ( this%time_start < 0.0_MK ) THEN
                 
                 PRINT *, "phsyics_check_parameters : ", &
                      "time_start does not match time_end !"
                 lcheck = .FALSE.
                 GOTO 9999
                 
              END IF
              
           END IF ! time_end >= 0.0
           
        END IF

        !----------------------------------------------------
        ! Check for fluid properties.
        !----------------------------------------------------

        !----------------------------------------------------
        ! Check for Brownian fluids, kt is given.
        !----------------------------------------------------
        
        IF ( Brownian .AND. &
             this%kt <= mcf_machine_zero) THEN
           PRINT *, "physics_check_parameters : ", &
                "Brownian fluid, kt should be given !"
           lcheck = .FALSE.
           GOTO 9999
        END IF

        IF ( this%rho_ref > this%rho + mcf_machine_zero ) THEN
           PRINT *, "physics_check_parameters : ", &
                "rho_ref should not be bigger than rho !"
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check for relax run, relax_* parameters are fine.
        !----------------------------------------------------
        
        IF ( relax_run ) THEN
           
           SELECT CASE(this%relax_type)

           CASE (1)
              
              IF ( this%step_relax <= 0 .AND. &
                   this%time_relax <= 0.0 ) THEN
                 PRINT *, "physics_check_parameters : ", &
                      "relax type 1, step_relax or time_relax should be positive !"
                 lcheck = .FALSE.
                 GOTO 9999
              END IF

              IF ( this%step_relax > 0 .AND. &
                   this%time_relax > 0.0_MK ) THEN
                 PRINT *, "phsyics_check_parameters : ", &
                      "relax type 1, only one of step and time should be non-negative !"
                 lcheck = .FALSE.
                 GOTO 9999
              END IF
        
           CASE (2)
              IF ( this%disorder_level <= 0.0_MK ) THEN
                 PRINT *, "physics_check_parameters : ", &
                      "relax run type 2, disorder should be positive !"
                 lcheck = .FALSE.
                 GOTO 9999
              END IF
              
           END SELECT
           
           IF ( this%kt_relax <= mcf_machine_zero) THEN
              PRINT *, "physics_check_parameters : ", &
                   "relax run, kt_relax should be given !"
              lcheck = .FALSE.
              GOTO 9999
           END IF
           
           IF ( this%c_relax <= mcf_machine_zero ) THEN
              PRINT *, "physics_check_parameters : ", &
                   "relax run, c_relax should be positive !"
              lcheck = .FALSE.
              GOTO 9999
           END IF
           
        END IF ! relax run
        
        !----------------------------------------------------
        ! Check viscoelastic parameters 
        !----------------------------------------------------

        IF ( .NOT. Newtonian ) THEN
        
           IF ( Brownian .AND. (.NOT. this%eigen_dynamics)) THEN
              PRINT *, "physics_check_parameters: ",&
                   "non-Newtonian Brownian should use eigen_dynamics !"
              lcheck = .FALSE.
              GOTO 9999
           END IF
           
           IF ( this%eigen_dynamics ) THEN
              
              IF ( this%evec_tolerance < mcf_machine_zero .OR. &
                   this%evec_tolerance > 1.0_MK-mcf_machine_zero ) THEN
                 PRINT *, "physics_check_parameters: ",&
                      "evec_tolerance is wrong ! "
                 lcheck = .FALSE.
                 GOTO 9999
              END IF
              
           END IF
           
        END IF

        !----------------------------------------------------
        ! Check if boundary definitions given resonable.
        !
        ! At most (2*num_dim -2) boundaries can be
        ! symmetry or wall boundaries.
        !----------------------------------------------------

        DO i = 1, dim
           
           IF ( this%bcdef(2*i-1) < 1 .OR. &
                this%bcdef(2*i-1) > 10 ) THEN
              PRINT *, "physics_check_parameters : ", &
                   "boundary condition  not defined !"
              lcheck = .FALSE.
              GOTO 9999
           END IF
           
           IF ( this%bcdef(2*i) < 1 .OR. &
                this%bcdef(2*i) > 10 ) THEN
              PRINT *, "physics_check_parameters : ", &
                   "boundary condition  not defined !"
              lcheck = .FALSE.
              GOTO 9999
           END IF
           
           IF ( this%bcdef(2*i-1) == ppm_param_bcdef_periodic .AND. &
                this%bcdef(2*i)  /= this%bcdef(2*i-1) ) THEN
              PRINT *, "physics_check_parameters : ", &
                   "periodic boundary must be pair-wise !"
              lcheck = .FALSE.
              GOTO 9999
           END IF

           IF ( this%bcdef(2*i) == ppm_param_bcdef_periodic .AND. &
                this%bcdef(2*i-1)  /= this%bcdef(2*i) ) THEN
              PRINT *, "physics_check_parameters : ", &
                   "periodic boundary must be pair-wise !"
              lcheck = .FALSE.
              GOTO 9999
           END IF

           IF ( this%bcdef(2*i-1) == ppm_param_bcdef_wall_solid .AND. &
                this%bcdef(2*i)  /= this%bcdef(2*i-1) ) THEN
              PRINT *, "physics_check_parameters : ", &
                   "solid wall boundary must be pair-wise !"
              lcheck = .FALSE.
              GOTO 9999
           END IF

           IF ( this%bcdef(2*i) == ppm_param_bcdef_wall_solid .AND. &
                this%bcdef(2*i-1)  /= this%bcdef(2*i) ) THEN
              PRINT *, "physics_check_parameters : ", &
                   "solid wall boundary must be pair-wise !"
              lcheck = .FALSE.
              GOTO 9999
           END IF

           IF ( this%bcdef(2*i-1) == ppm_param_bcdef_LE .AND. &
                this%bcdef(2*i)  /= this%bcdef(2*i-1) ) THEN
              PRINT *, "physics_check_parameters : ", &
                   "Lees-Edwards boundary must be pair-wise !"
              lcheck = .FALSE.
              GOTO 9999
           END IF
           
           IF ( this%bcdef(2*i) == ppm_param_bcdef_LE .AND. &
                this%bcdef(2*i-1)  /= this%bcdef(2*i) ) THEN
              PRINT *, "physics_check_parameters : ", &
                   "Lees-Edwards boundary must be pair-wise !"
              lcheck = .FALSE.
              GOTO 9999
           END IF
           
        END DO ! i = 1, num_dim
        
        
        IF ( symmetry .AND. num_sym > 0 ) THEN
           
           PRINT *, "physics_check_parameters : ", &
                "Symmetry inter-process, no symmetry boundary !"
           lcheck = .FALSE.
           GOTO 9999
           
        END IF
        
        IF ( symmetry .AND. num_wall_sym > 0 ) THEN
           
           PRINT *, "phsyics_check_parameters : ", &
                "Symmetry inter-process, no wall_sym boundary !"
                lcheck = .FALSE.
           GOTO 9999       
           
        END IF
        
        IF ( num_sym > 2*this%num_dim ) THEN
           
           PRINT *, "phsyics_check_parameters : ", &
                "At most 2*num_dim - 2 boundaries can be symmetry"
           lcheck = .FALSE.
           GOTO 9999
           
        END IF
        
        
        IF ( num_wall_sym > 0 .AND. &
             num_wall_solid > 0 ) THEN
           
           PRINT *, "phsyics_check_parameters : ", &
                "Either wall_sym by PPM or, wall_solid by MCF !"
           lcheck = .FALSE.
           GOTO 9999         
           
        END IF
        
        
#if 0
        IF ( num_wall > 2*this%num_dim - 2 ) THEN
           
           PRINT *, "phsyics_check_parameters : ", &
                "At most 2*num_dim - 2 boundaries can be wall"
           lcheck = .FALSE.
           GOTO 9999         
           
        END IF
#endif   
        
        !----------------------------------------------------
        ! If number of species bigger than 1,
        ! we have to have at least one colloid.
        !----------------------------------------------------
        
        IF ( this%num_species > 1  .AND. &
             this%num_colloid < 1) THEN
           PRINT *, "physics_check_parameters : ", &
                "At least one colloid in 2 species case !"
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! For two species with colloids,
        ! we check if colloids parameters reasonable.
        !----------------------------------------------------
        
        
        IF ( this%num_species > 1  .AND. &
             this%num_colloid > 0) THEN
           
           !-------------------------------------------------
           ! Check further colloids parameters.
           !-------------------------------------------------
           
           lcheck = &
                colloid_check_parameters(this%colloids,stat_info_sub)
           
        END IF
        
9999    CONTINUE
        
        physics_check_parameters = lcheck
        
        RETURN
        
      END FUNCTION physics_check_parameters
