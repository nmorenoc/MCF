      SUBROUTINE io_write_condition_check(this,&
           step_current,time_current,&
           write_particles, write_conformation,&
           write_colloid, &
           write_statistic, write_boundary,&
           write_restart_physics,&
           write_restart_particles,&
           write_restart_conformation, stat_info)
        
        !----------------------------------------------------
        !  Subroutine   :  io_write_condition_check
        !----------------------------------------------------
        !
        !  Purpose      :  Check if write or not at current
        !                  step/time/wall_time
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : 
        !                 V0.1 20.01.2009, original version.
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
        
        !----------------------------------------------------
        ! Arguments.
        !----------------------------------------------------
        
        TYPE(IO), INTENT(INOUT)         :: this
        INTEGER, INTENT(IN)             :: step_current
        REAL(MK), INTENT(IN)            :: time_current
        LOGICAL, INTENT(OUT),OPTIONAL   :: write_particles
        LOGICAL, INTENT(OUT),OPTIONAL   :: write_conformation
        LOGICAL, INTENT(OUT),OPTIONAL   :: write_colloid
        LOGICAL, INTENT(OUT),OPTIONAL   :: write_statistic
        LOGICAL, INTENT(OUT),OPTIONAL   :: write_boundary
        LOGICAL, INTENT(OUT),OPTIONAL   :: write_restart_physics
        LOGICAL, INTENT(OUT),OPTIONAL   :: write_restart_particles
        LOGICAL, INTENT(OUT),OPTIONAL   :: write_restart_conformation        
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: step
        REAL(MK)                        :: time
        REAL(MK)                        :: wall_time
        LOGICAL                         :: Newtonian
        INTEGER                         :: num_colloid
        INTEGER                         :: write_restart
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: num_shear
         
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        NULLIFY(tboundary)
        
        !----------------------------------------------------
        ! Step since simulation starts at current run.
        ! Time since simulation starts at current run.
        !----------------------------------------------------
        
        step = step_current - this%step_start
        time = time_current - this%time_start
        
        !----------------------------------------------------
        ! Get control and physical parameters.
        !----------------------------------------------------
        
        Newtonian     = &
             control_get_Newtonian(this%ctrl,stat_info_sub)
        num_colloid   = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        Write_restart = this%write_restart
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        num_shear     = &
             boundary_get_num_shear(tboundary,stat_info_sub)
        
        
        !----------------------------------------------------
        ! Checking writting for output according to different
        ! time measurement,
        ! 1: simulation step
        ! 2: simulation time
        !----------------------------------------------------

        SELECT CASE (this%write_output)
           
           !-------------------------------------------------
           ! Writting freq is based on simulation step.
           !-------------------------------------------------
           
        CASE (0)
           
           this%write_particles    = .FALSE.
           this%write_conformation = .FALSE.
           this%write_colloid      = .FALSE.
           this%write_statistic    = .FALSE.
           this%write_boundary     = .FALSE.
           
        CASE (1)

           IF(MOD(step,this%output_particles_freq_step) == 0 ) THEN
              this%write_particles = .TRUE.
           ELSE
              this%write_particles = .FALSE.
           END IF
           
           IF( (.NOT. Newtonian ).AND. &
                MOD(step,this%output_conformation_freq_step) == 0 ) THEN
              this%write_conformation = .TRUE.
           ELSE
              this%write_conformation = .FALSE.
           END IF
           
           IF( num_colloid > 0 .AND. &
                MOD(step,this%colloid_freq_step) ==0 ) THEN
              this%write_colloid = .TRUE.
           ELSE
              this%write_colloid = .FALSE.
           END IF

           IF(MOD(step, this%statistic_freq_step) ==0 ) THEN
              this%write_statistic = .TRUE.
           ELSE
              this%write_statistic = .FALSE.
           END IF
           
           IF( ( num_shear > 0 )  .AND.&
                MOD(step, this%boundary_freq_step) == 0 ) THEN
              this%write_boundary = .TRUE.
           ELSE
              this%write_boundary = .FALSE.
           END IF
           
           !-------------------------------------------------
           ! Writting freq is based on simulation time.
           !-------------------------------------------------
           
        CASE (2)
           
           IF( time >= this%output_particles_freq_time * &
                this%output_particles_freq_time_num ) THEN
              this%write_particles = .TRUE.
              
              this%output_particles_freq_time_num = &
                   this%output_particles_freq_time_num + 1
           ELSE
              this%write_particles = .FALSE.
           END IF
           
           IF( (.NOT. Newtonian ).AND. &
                time >= this%output_conformation_freq_time * &
                this%output_conformation_freq_time_num ) THEN
              
              this%write_conformation = .TRUE.
              this%output_conformation_freq_time_num = &
                   this%output_conformation_freq_time_num + 1
           ELSE
              
              this%write_conformation = .FALSE.

           END IF
           
           IF( num_colloid > 0 .AND. &
                time >=this%colloid_freq_time * &
                this%colloid_freq_time_num ) THEN
              
              this%write_colloid = .TRUE.
              this%colloid_freq_time_num = &
                   this%colloid_freq_time_num + 1
              
           ELSE
              
              this%write_colloid = .FALSE.
              
           END IF
           
           IF( time >= this%statistic_freq_time * &
                this%statistic_freq_time_num ) THEN
              
              this%write_statistic = .TRUE.
              this%statistic_freq_time_num = &
                   this%statistic_freq_time_num + 1
              
           ELSE
              
              this%write_statistic = .FALSE.
              
           END IF
           
           IF( ( num_shear > 0 )  .AND.&
                time >= this%boundary_freq_time * &
                this%boundary_freq_time_num ) THEN
              
              this%write_boundary = .TRUE.
              this%boundary_freq_time_num = &
                   this%boundary_freq_time_num + 1
              
           ELSE
              
              this%write_boundary = .FALSE.
              
           END IF
           
        CASE DEFAULT
           
           PRINT *,__FILE__, __LINE__, &
                "no such way of writing output files!"
           stat_info=-1
           GOTO 9999
           
        END SELECT ! write_output
        
        !----------------------------------------------------
        ! Checking writting for restart according to different
        ! time measurement,
        ! 1: simulation step
        ! 2: simulation time
        ! 3: wall time
        !----------------------------------------------------

        
        SELECT CASE ( this%write_restart ) 
           
        CASE (0)
           this%write_restart_physics      = .FALSE.
           this%write_restart_particles    = .FALSE.
           this%write_restart_conformation = .FALSE.
                 
        CASE (1)
           
           IF( step > 0 .AND. &
                MOD(step,this%restart_freq_step) == 0 ) THEN
              
              this%write_restart_physics   = .TRUE.           
              this%write_restart_particles = .TRUE.
              
              IF ( .NOT. Newtonian ) THEN
                 this%write_restart_conformation = .TRUE.
              ELSE
                 this%write_restart_conformation = .FALSE.
              END IF
              
           ELSE
              
              this%write_restart_physics      = .FALSE. 
              this%write_restart_particles    = .FALSE.
              this%write_restart_conformation = .FALSE.
              
           END IF
           
        CASE (2)
           
           IF( time >= this%restart_freq_time * &
                this%restart_freq_time_num  ) THEN
              
              this%restart_freq_time_num = &
                   this%restart_freq_time_num + 1
              
              this%write_restart_physics   = .TRUE.           
              this%write_restart_particles = .TRUE.
              
              IF ( .NOT. Newtonian ) THEN
                 this%write_restart_conformation = .TRUE.
              ELSE
                 this%write_restart_conformation = .FALSE.
              END IF
              
           ELSE
              
              this%write_restart_physics      = .FALSE.
              this%write_restart_particles    = .FALSE.
              this%write_restart_conformation = .FALSE.
              
           END IF
           
        CASE (3)

           wall_time = &
                control_get_elapsed_wall_time(this%ctrl,stat_info_sub)
           
           IF( wall_time >= this%restart_freq_time_wall * &
                this%restart_freq_time_num  ) THEN
              
              this%restart_freq_time_num = &
                   this%restart_freq_time_num + 1
              
              this%write_restart_physics   = .TRUE.           
              this%write_restart_particles = .TRUE.
              
              IF ( .NOT. Newtonian ) THEN
                 this%write_restart_conformation = .TRUE.
              ELSE
                 this%write_restart_conformation = .FALSE.
              END IF
              
           ELSE
              
              this%write_restart_physics      = .FALSE.
              this%write_restart_particles    = .FALSE.
              this%write_restart_conformation = .FALSE.
              
           END IF
           
           
        CASE DEFAULT
           
           PRINT *,__FILE__, __LINE__, &
                "no such way of writing restart files!"
           stat_info=-1
           GOTO 9999

        END SELECT ! write_restart
         
        
        IF (PRESENT(write_particles) ) THEN
           write_particles = this%write_particles
        END IF
        
        IF (PRESENT(write_conformation) ) THEN
           write_conformation = this%write_conformation
        END IF
        
        IF (PRESENT(write_colloid) ) THEN
           write_colloid = this%write_colloid
        END IF
      
        IF (PRESENT(write_statistic) ) THEN
           write_statistic = this%write_statistic
        END IF
        
        IF (PRESENT(write_boundary) ) THEN
           write_boundary = this%write_boundary
        END IF
        
        
        IF (PRESENT(write_restart_physics) ) THEN
           write_restart_physics = this%write_restart_physics
        END IF
        
        IF (PRESENT(write_restart_particles) ) THEN
           write_restart_particles = this%write_restart_particles
        END IF
        
        IF (PRESENT(write_restart_conformation) ) THEN
           write_restart_conformation = this%write_restart_conformation
        END IF
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE io_write_condition_check
      
      
      SUBROUTINE io_write_condition_set(this,&
           write_particles, write_conformation,&
           write_colloid, &
           write_statistic, write_boundary,&          
           write_restart_physics,&
           write_restart_particles,&
           write_restart_conformation, stat_info)
        
        !----------------------------------------------------
        !  Subroutine   :  io_write_condition_set
        !----------------------------------------------------
        !
        !  Purpose      :  set write condition
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : 
        !                 V0.1 13.10.2009, original version.
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
        
        !----------------------------------------------------
        ! Arguments
        !----------------------------------------------------
        
        TYPE(IO), INTENT(INOUT)         :: this
        LOGICAL, INTENT(IN),OPTIONAL    :: write_particles
        LOGICAL, INTENT(IN),OPTIONAL    :: write_conformation
        LOGICAL, INTENT(IN),OPTIONAL    :: write_colloid
        LOGICAL, INTENT(IN),OPTIONAL    :: write_statistic
        LOGICAL, INTENT(IN),OPTIONAL    :: write_boundary        
        LOGICAL, INTENT(IN),OPTIONAL    :: write_restart_physics
        LOGICAL, INTENT(IN),OPTIONAL    :: write_restart_particles
        LOGICAL, INTENT(IN),OPTIONAL    :: write_restart_conformation        
        INTEGER, INTENT(INOUT)          :: stat_info
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: Newtonian
        INTEGER                         :: num_colloid
        INTEGER                         :: write_restart
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: num_shear
        
        stat_info = 0
        NULLIFY(tboundary)
        
        Newtonian     = &
             control_get_Newtonian(this%ctrl,stat_info_sub)
        num_colloid   = &
             physics_get_num_colloid(this%phys,stat_info_sub)        
        Write_restart = this%write_restart
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        num_shear     = &
             boundary_get_num_shear(tboundary,stat_info_sub)
        
        IF (PRESENT(write_particles) ) THEN
           this%write_particles = write_particles
        END IF
        
        IF (PRESENT(write_conformation) .AND. &
             (.NOT. Newtonian)) THEN
           this%write_conformation = write_conformation
        END IF
        
        IF ( num_colloid >0 .AND. &
             PRESENT(write_colloid) ) THEN
           this%write_colloid =  write_colloid
        END IF
        
        IF (PRESENT(write_statistic) ) THEN
           this%write_statistic = write_statistic
        END IF
        
        IF ( ( num_shear > 0 )  .AND. &
             PRESENT(write_boundary)) THEN
           this%write_boundary =  write_boundary
        END IF
        
        IF ( write_restart > 0 .AND. &
             PRESENT(write_restart_physics) ) THEN
           this%write_restart_physics = write_restart_physics
        END IF
        
        IF ( write_restart > 0 .AND. &
             PRESENT(write_restart_particles) ) THEN
           this%write_restart_particles = write_restart_particles
        END IF
        
        IF ( write_restart > 0.AND. &
             (.NOT. Newtonian ) .AND. &
             PRESENT(write_restart_conformation) ) THEN
           this%write_restart_conformation = write_restart_conformation
        END IF
        
        RETURN
        
      END SUBROUTINE io_write_condition_set
      
      
