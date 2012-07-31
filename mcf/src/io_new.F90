      SUBROUTINE io_init_default(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_init_default
        !----------------------------------------------------
        !
        ! Purpose     : Default construtor of IO Class.
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

        !----------------------------------------------------
        ! Arguments
        !----------------------------------------------------

        TYPE(IO),INTENT(OUT)            :: this
        INTEGER,INTENT(OUT)             :: stat_info
        
        INTEGER                         :: stat_info_sub
        
        stat_info     = 0
        stat_info_sub = 0
        
        this%ctrl_file = TRIM("ctrl.mcf")
        this%ctrl_unit = 20
        
        this%physics_config_file = TRIM("physics_config.mcf")
        this%physics_config_unit = 21
        
        this%io_config_file = TRIM("io_config.mcf")
        this%io_config_unit = 22
        
        this%write_output = 1
        
        this%read_particles_file  = TRIM("mcf_restart_particles.dat")
        this%read_particles_unit  = 23
        this%read_particles_fmt   =  TRIM("FORMATTED")

        this%read_conformation_file  = TRIM("mcf_restart_conformation.dat")
        this%read_conformation_unit  = 24
        this%read_conformation_fmt   =  TRIM("FORMATTED")
        
        this%output_particles_relax_file  = TRIM("mcf_init_particles")
        this%output_particles_relax_unit  = 30
        this%output_particles_relax_fmt   = TRIM("FORMATTED")
        this%output_particles_relax_freq_step  = 100

        this%output_particles_file  = TRIM("mcf_particles")
        this%output_particles_unit  = 31
        this%output_particles_fmt   = TRIM("FORMATTED")
        this%output_particles_freq_step  = 100
        this%output_particles_freq_time  = 0.1_MK
        this%output_particles_freq_time_num = 0

        this%output_conformation_file  = TRIM("mcf_conformation")
        this%output_conformation_unit  = 32
        this%output_conformation_fmt   = TRIM("FORMATTED")
        this%output_conformation_freq_step  = 100
        this%output_conformation_freq_time  = 0.1_MK
        this%output_conformation_freq_time_num = 0

        this%colloid_file  = TRIM("mcf_colloid")
        this%colloid_unit  = 100
        this%colloid_fmt   = TRIM("FORMATTED")
        this%colloid_freq_step  = 100
        this%colloid_freq_time  = 0.1_MK
        this%colloid_freq_time_num = 0
        
        this%statistic_relax_file       = TRIM("mcf_init_statistic.dat")
        this%statistic_relax_unit       = 40
        this%statistic_relax_fmt        = TRIM("FORMATTED")
        this%statistic_relax_freq_step  = 1

        this%statistic_file  = TRIM("mcf_statistic.dat")
        this%statistic_unit  = 41
        this%statistic_fmt   = TRIM("FORMATTED")
        this%statistic_freq_step  = 100
        this%statistic_freq_time  = 0.1_MK
        this%statistic_freq_time_num = 0

        this%boundary_file  = TRIM("mcf_boundary.dat")
        this%boundary_unit  = 50
        this%boundary_fmt   = TRIM("FORMATTED")
        this%boundary_freq_step  = 100
        this%boundary_freq_time  = 0.1_MK
        this%boundary_freq_time_num = 0
        
        this%write_restart = 0
        
        this%restart_particles_relax_file  = TRIM("mcf_init_restart_particles")
        this%restart_particles_relax_unit  = 90
        this%restart_particles_relax_fmt   =  TRIM("FORMATTED")
        
        this%restart_physics_file    = TRIM("mcf_restart_physics")
        this%restart_physics_unit    = 91
        this%restart_physics_fmt     =  TRIM("FORMATTED")
        
        this%restart_particles_file  = TRIM("mcf_restart_particles")
        this%restart_particles_unit  = 92
        this%restart_particles_fmt   =  TRIM("FORMATTED")
        
        this%restart_conformation_file  = TRIM("mcf_restart_conformation")
        this%restart_conformation_unit  = 93
        this%restart_conformation_fmt   =  TRIM("FORMATTED")
        
        this%restart_freq_step       = 100
        this%restart_freq_time       = 0.1_MK
        this%restart_freq_time_wall  = 48.0_MK
        this%restart_freq_time_num   = 1

        this%write_particles    = .FALSE.
        this%write_conformation = .FALSE.
        this%write_colloid      = .FALSE.
        this%write_statistic    = .FALSE.
        this%write_boundary     = .FALSE.        
        this%write_restart_physics      = .FALSE.
        this%write_restart_particles    = .FALSE.
        this%write_restart_conformation = .FALSE.

        CALL tool_new(this%io_tool,stat_info_sub)
        
        
        RETURN          
        
      END SUBROUTINE io_init_default
      
      
      SUBROUTINE io_init(this,d_ctrl,d_ctrl_file,d_phys,stat_info)
        
        TYPE(IO), INTENT(OUT)                  :: this
        TYPE(Control), INTENT(INOUT),TARGET    :: d_ctrl
        CHARACTER(LEN=*), INTENT(IN)           :: d_ctrl_file
        TYPE(Physics), INTENT(INOUT),TARGET    :: d_phys
        INTEGER,INTENT(OUT)                    :: stat_info
        
        INTEGER                                :: stat_info_sub
        
        stat_info     = 0
        stat_info_sub = 0
        
        !------------------------------------------
        ! The name of control file is by default 
        ! "ctrl.mcf", if not given as parameter.
        ! The other two files which are physics 
        ! config and io config files, named in
        ! in "ctrl.mcf".
        !--------------------------------------------------------------

        !--------------------------------------------------------------
        ! The files' Units of file handling are all pre-defined, fixed.
        !--------------------------------------------------------------

        NULLIFY(this%ctrl)
        this%ctrl => d_ctrl
        
        NULLIFY(this%phys)
        this%phys => d_phys
     
        IF ( LEN_TRIM(d_ctrl_file) < 0 ) THEN
           PRINT *, "io_init : ", &
                "The given ctrl file name is empty !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%ctrl_file = TRIM(ADJUSTL(d_ctrl_file))
        
        this%ctrl_unit = 20
        
        this%physics_config_file = TRIM("physics_config.mcf")
        this%physics_config_unit = 21
        
        this%io_config_file = TRIM("io_config.mcf")
        this%io_config_unit = 22      
        
        this%read_particles_file  = TRIM("mcf_restart_particles.dat")
        this%read_particles_unit  = 23
        this%read_particles_fmt   =  TRIM("FORMATTED")
        
        this%read_conformation_file  = TRIM("mcf_restart_conformation.dat")
        this%read_conformation_unit  = 24
        this%read_conformation_fmt   =  TRIM("FORMATTED")
        
        this%output_particles_relax_file  = TRIM("mcf_init_particles")
        this%output_particles_relax_unit  = 30
        this%output_particles_relax_fmt   =  TRIM("FORMATTED")
        this%output_particles_relax_freq_step  = 100
        
        this%output_particles_file  = TRIM("mcf_particles")
        this%output_particles_unit  = 31
        this%output_particles_fmt   =  TRIM("FORMATTED")
        this%output_particles_freq_step  = 100
        this%output_particles_freq_time  = 0.1_MK
        this%output_particles_freq_time_num = 0

        this%output_conformation_file  = TRIM("mcf_conformation")
        this%output_conformation_unit  = 32
        this%output_conformation_fmt   =  TRIM("FORMATTED")
        this%output_conformation_freq_step  = 100
        this%output_conformation_freq_time  = 0.1_MK
        this%output_conformation_freq_time_num = 0

        this%colloid_file  = TRIM("mcf_colloid")
        this%colloid_unit  = 100
        this%colloid_fmt   =  TRIM("FORMATTED")
        this%colloid_freq_step  = 1
        this%colloid_freq_time  = 0.1_MK
        this%colloid_freq_time_num = 0

        this%statistic_relax_file  =  TRIM("mcf_init_statistic.dat")
        this%statistic_relax_unit  = 40
        this%statistic_relax_fmt   =  TRIM("FORMATTED")
        this%statistic_relax_freq_step  = 1
        
        this%statistic_file  = TRIM("mcf_statistic.dat")
        this%statistic_unit  = 41
        this%statistic_fmt   =  TRIM("FORMATTED")
        this%statistic_freq_step  = 1
        this%statistic_freq_time  = 0.1_MK
        this%statistic_freq_time_num = 0

        this%boundary_file  =  TRIM("mcf_boundary.dat")
        this%boundary_unit  = 50
        this%boundary_fmt   =  TRIM("FORMATTED")
        this%boundary_freq_step  = 1
        this%boundary_freq_time  = 0.1_MK
        this%boundary_freq_time_num = 0

        this%restart_particles_relax_file  = TRIM("mcf_init_restart_particles")
        this%restart_particles_relax_unit  = 90
        this%restart_particles_relax_fmt   =  TRIM("FORMATTED")

        this%restart_physics_file    = TRIM("mcf_restart_physics")
        this%restart_physics_unit    = 91
        this%restart_physics_fmt     =  TRIM("FORMATTED")
        
        this%restart_particles_file  = TRIM("mcf_restart_particles")
        this%restart_particles_unit  = 92
        this%restart_particles_fmt   =  TRIM("FORMATTED")
        
        this%restart_conformation_file  = TRIM("mcf_restart_conformation")
        this%restart_conformation_unit  = 93
        this%restart_conformation_fmt   =  TRIM("FORMATTED")
        
        this%restart_freq_step  = 100
        this%restart_freq_time  = 0.1_MK
        this%restart_freq_time_wall = 48.0_MK
        this%restart_freq_time_num   = 1

        CALL tool_new(this%io_tool,stat_info_sub)        
        
        CALL io_read_ctrl(this,d_ctrl,stat_info_sub)
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "io_init : ", &
                "Reading control file has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        CALL io_read_physics_config(this,d_phys,stat_info_sub)
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "io_init : ", &
                "Reading physics config file has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        CALL io_read_io_config(this,stat_info_sub)
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "io_init : ", &
                "Reading io config file has problem !"
           stat_info = -1
           GOTO 9999
        END IF

        this%write_output  = &
             control_get_write_output(this%ctrl,stat_info_sub)
        this%write_restart = &
             control_get_write_restart(this%ctrl,stat_info_sub)
        this%write_particles    = .FALSE.
        this%write_conformation = .FALSE.
        this%write_colloid      = .FALSE.
        this%write_statistic    = .FALSE.
        this%write_boundary     = .FALSE.        
        this%write_restart_physics      = .FALSE.
        this%write_restart_particles    = .FALSE.
        this%write_restart_conformation = .FALSE.
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE io_init
      
      
      SUBROUTINE io_display_parameters(this,stat_info)
        
        
        TYPE(IO),INTENT(IN)           :: this
        INTEGER,INTENT(OUT)           :: stat_info
        
        LOGICAL                       :: read_external
        LOGICAL                       :: relax_run
        LOGICAL                       :: Newtonian
        
        INTEGER                       :: num_species
        TYPE(Boundary), POINTER       :: tboundary
        INTEGER                       :: num_shear
        INTEGER                       :: num_colloid

        INTEGER                       :: stat_info_sub
        
        stat_info     = 0
        stat_info_sub = 0
        NULLIFY(tboundary)
        
        read_external = &
             control_get_read_external(this%ctrl,stat_info_sub)
        relax_run     = &
             control_get_relax_run(this%ctrl,stat_info_sub)  
        Newtonian     = &
             control_get_Newtonian(this%ctrl,stat_info_sub)
        num_species   = &
             physics_get_num_species(this%phys,stat_info_sub)
        
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        
        num_shear    = &
             boundary_get_num_shear(tboundary,stat_info_sub)
        
        num_colloid  = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        PRINT *, '------------------Start------------------'
        PRINT *, '     IO parameters'
        PRINT *, '-----------------------------------------'
        
        PRINT *, "ctrl_file                   : ", TRIM(this%ctrl_file)
        PRINT *, "ctrl_unit                   : ", this%ctrl_unit
        
        PRINT *, "physics_config_file         : ", TRIM(this%physics_config_file)
        PRINT *, "physics_config_unit         : ", this%physics_config_unit
        
        PRINT *, "io_config_file              : ", TRIM(this%io_config_file)
        PRINT *, "io_config_unit              : ", this%io_config_unit
        
        IF( read_external ) THEN
           
           PRINT *, "read_particle_file          : ", TRIM(this%read_particles_file)
           PRINT *, "read_particles_unit         : ", this%read_particles_unit
           PRINT *, "read_particles_format       : ", TRIM(this%read_particles_fmt)
           
           IF ( .NOT. Newtonian ) THEN
              
              PRINT *, "read_conformation_file      : ", TRIM(this%read_conformation_file)
              PRINT *, "read_conformation_unit      : ", this%read_conformation_unit
              PRINT *, "read_conformation_format    : ", TRIM(this%read_conformation_fmt)
              
           END IF
           
        END IF
        
        PRINT *, "write output                : ", this%write_output
        IF ( this%write_output > 0 .AND. &
             this%write_output < 3 ) THEN
           
           IF( relax_run ) THEN
              PRINT *, "output_particles_relax_file      : ", TRIM(this%output_particles_relax_file), &
                   "****.out"
              PRINT *, "output_particles_relax_unit      : ", this%output_particles_relax_unit
              PRINT *, "output_particles_relax_format    : ", TRIM(this%output_particles_relax_fmt)
              PRINT *, "output_particles_relax_freq_step : ", this%output_particles_relax_freq_step
           END IF
           
           PRINT *, "output_particles_file      : ", TRIM(this%output_particles_file), &
                "****.out"
           PRINT *, "output_particles_unit      : ", this%output_particles_unit
           PRINT *, "output_particles_format    : ", TRIM(this%output_particles_fmt)
           
           SELECT CASE ( this%write_output )
           CASE (1)
              PRINT *, "output_particles_freq_step : ", this%output_particles_freq_step
           CASE (2)
              PRINT *, "output_particles_freq_time : ", this%output_particles_freq_time
           END SELECT
           
           IF ( .NOT. Newtonian ) THEN
              PRINT *, "output_conformation_file      : ", TRIM(this%output_conformation_file), &
                   "****.out"
              PRINT *, "output_conformation_unit      : ", this%output_conformation_unit
              PRINT *, "output_conformation_format    : ", TRIM(this%output_conformation_fmt)
              SELECT CASE ( this%write_output )
              CASE (1)
                 PRINT *, "output_conformation_freq_step : ", this%output_conformation_freq_step
              CASE (2)
                 PRINT *, "output_conformation_freq_time : ", this%output_conformation_freq_time
              END SELECT
           END IF
           
           IF ( num_species > 1 .AND. num_colloid > 0 )THEN
              PRINT *, "colloid_file      : ", TRIM(this%colloid_file), &
                   "**.dat"
              PRINT *, "colloid_unit      : ", this%colloid_unit
              PRINT *, "colloid_format    : ", TRIM(this%colloid_fmt)
              SELECT CASE ( this%write_output )
              CASE (1)
                 PRINT *, "colloid_freq_step : ", this%colloid_freq_step
              CASE (2)
                 PRINT *, "colloid_freq_time : ", this%colloid_freq_time
              END SELECT
           END IF
           
           IF( relax_run ) THEN
              PRINT *, "statistic_relax_file      : ", TRIM(this%statistic_relax_file)
              PRINT *, "statistic_relax_unit      : ", this%statistic_relax_unit
              PRINT *, "statistic_relax_format    : ", TRIM(this%statistic_relax_fmt)
              PRINT *, "statistic_relax_freq_step : ", this%statistic_relax_freq_step
           END IF
           
           PRINT *, "statistic_file      : ", TRIM(this%statistic_file)
           PRINT *, "statistic_unit      : ", this%statistic_unit
           PRINT *, "statistic_format    : ", TRIM(this%statistic_fmt)
           SELECT CASE ( this%write_output )
           CASE (1)
              PRINT *, "statistic_freq_step : ", this%statistic_freq_step
           CASE (2)
              PRINT *, "statistic_freq_time : ", this%statistic_freq_time
           END SELECT
           
           IF ( num_shear > 0 )THEN
              PRINT *, "boundary_file      : ", TRIM(this%boundary_file)
              PRINT *, "boundary_unit      : ", this%boundary_unit
              PRINT *, "boundary_format    : ", TRIM(this%boundary_fmt)
              SELECT CASE ( this%write_output )
              CASE (1)
                 PRINT *, "boundary_freq_step : ", this%boundary_freq_step
              CASE (2)
                 PRINT *, "boundary_freq_time : ", this%boundary_freq_time
              END SELECT
           END IF
           
           IF( relax_run ) THEN
              IF( LEN(TRIM(this%restart_particles_relax_file)) > 0 ) THEN
                 PRINT *, "restart_particles_relax_file   : ", TRIM(this%restart_particles_relax_file),&
                      "****.dat"
              ELSE
                 PRINT *, "restart_particles_relax_file   : "
              END IF
              PRINT *, "restart_particles_relax_unit   : ", this%restart_particles_relax_unit
              PRINT *, "restart_particles_relax_format : ", TRIM(this%restart_particles_relax_fmt)
              
           END IF
           
        END IF ! write_output
        
        
        IF ( this%write_restart > 0 ) THEN
           
           IF(LEN(TRIM(this%restart_physics_file)) > 0) THEN
              PRINT *, "restart_physics_file   : ", TRIM(this%restart_physics_file),&
                   "****.dat"
           ELSE
              PRINT *, "restart_physics_file   : "
           END IF
           
           PRINT *, "restart_physics_unit   : ", this%restart_physics_unit
           PRINT *, "restart_physics_format : ", TRIM(this%restart_physics_fmt)
           
           
           IF(LEN(TRIM(this%restart_particles_file)) > 0) THEN
              PRINT *, "restart_particles_file   : ", TRIM(this%restart_particles_file),&
                   "****.dat"
           ELSE
              PRINT *, "restart_particles_file   : "
           END IF
              
           PRINT *, "restart_particles_unit   : ", this%restart_particles_unit
           PRINT *, "restart_particles_format : ", TRIM(this%restart_particles_fmt)
           
           IF ( .NOT. Newtonian ) THEN
              
              IF(LEN(TRIM(this%restart_conformation_file)) > 0) THEN
                 PRINT *, "restart_conformation_file   : ", TRIM(this%restart_conformation_file),&
                      "****.dat"
              ELSE
                 PRINT *, "restart_conformation_file   : "
              END IF
              
              PRINT *, "restart_conformation_unit   : ", this%restart_conformation_unit
              PRINT *, "restart_conformation_format : ", TRIM(this%restart_conformation_fmt)
              
           END IF
           
           SELECT CASE ( this%write_restart ) 
           CASE (1)
              PRINT *, "restart_freq_step        : ", this%restart_freq_step
           CASE (2)
              PRINT *, "restart_freq_time        : ", this%restart_freq_time
           CASE (3)
              PRINT *, "restart_freq_time_wall   : ", &
                   this%restart_freq_time_wall, " hours"
           CASE DEFAULT
              PRINT *, "no such way of writting restart files!"
           END SELECT
           
        END IF ! write_restart
        
        
        PRINT *, '-------------------End-------------------'
        
        
        RETURN
        
      END SUBROUTINE io_display_parameters
      
      
