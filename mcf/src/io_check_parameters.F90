      LOGICAL FUNCTION io_check_parameters(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_check_parameters
        !----------------------------------------------------
        !
        ! Purpose     :  Check io parameters to see if
        !                they are reasonable.
        !
        !  Revision   : V0.1 01.04.2009
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
        ! Argument
        !----------------------------------------------------

        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local parameters.
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: relax_run
        LOGICAL                         :: Newtonian
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: num_shear
        INTEGER                         :: num_colloid        
      

        !----------------------------------------------------
        ! Initialization.
        !----------------------------------------------------
        
        stat_info = 0
        io_check_parameters = .TRUE.
        
        relax_run = control_get_relax_run(this%ctrl,stat_info_sub)
        Newtonian = control_get_Newtonian(this%ctrl,stat_info_sub)
        
        NULLIFY(tboundary)
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        num_shear   = boundary_get_num_shear(tboundary,stat_info_sub)
        num_colloid = physics_get_num_colloid(this%phys,stat_info_sub)
        
        
        SELECT CASE(this%write_output)
           
        CASE (1)
           
           IF ( relax_run .AND. &
                this%output_particles_relax_freq_step <= 0 ) THEN
              PRINT *, "io_check_parameters : ", &
                   "output_particles_relax_freq_step <=0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF
           
           IF ( this%output_particles_freq_step <= 0 ) THEN
              PRINT *, "io_check_parameters : ", &
                   "output_particles_freq_step <=0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF
           
           IF ( .NOT. Newtonian .AND. &
                this%output_conformation_freq_step <= 0 ) THEN
              PRINT *, "io_check_parameters : ", &
                   "output_conformation_freq_step <=0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF

           IF ( num_colloid > 0 .AND. &
                this%colloid_freq_step <= 0 ) THEN
              PRINT *, "io_check_parameters : ", &
                   "colloid_freq_step <=0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF

           IF ( relax_run .AND. &
                this%statistic_relax_freq_step <= 0 ) THEN
              PRINT *, "io_check_parameters : ", &
                   "statistic_relax_freq_step <=0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF

           IF ( this%statistic_freq_step <= 0 ) THEN
              PRINT *, "io_check_parameters : ", &
                   "statistic_freq_step <=0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF
           
           IF ( num_shear > 0 .AND. &
                this%boundary_freq_step <= 0 ) THEN
              PRINT *, "io_check_parameters : ", &
                   "boundary_freq_step <=0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF
           
           
        CASE (2)
           
           IF ( this%output_particles_freq_time <= 0.0_MK ) THEN
              PRINT *, "io_check_parameters : ", &
                   "output_particles_freq_time <=0.0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF
           
           IF ( .NOT. Newtonian .AND. &
                this%output_conformation_freq_time <= 0.0_MK ) THEN
              PRINT *, "io_check_parameters : ", &
                   "output_conformation_freq_time <=0.0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF

           IF ( num_colloid > 0 .AND. &
                this%colloid_freq_time <= 0.0_MK ) THEN
              PRINT *, "io_check_parameters : ", &
                   "colloid_freq_time <=0.0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF
           
           IF ( this%statistic_freq_time <= 0.0_MK ) THEN
              PRINT *, "io_check_parameters : ", &
                   "statistic_freq_time <=0.0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF
           
           IF ( num_shear > 0 .AND. &
                this%boundary_freq_time <= 0.0_MK ) THEN
              PRINT *, "io_check_parameters : ", &
                   "boundary_freq_time <=0.0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF
           
        END SELECT ! wirte_output
        
        
        SELECT CASE (this%write_restart)
           
        CASE (1)
           
           IF ( this%restart_freq_step <= 0 ) THEN
              PRINT *, "io_check_parameters : ", &
                   "restart_freq_step <=0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF
           
        CASE (2)
           
           IF ( this%restart_freq_time <= 0.0_MK ) THEN
              PRINT *, "io_check_parameters : ", &
                   "restart_freq_time <=0.0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF
           
        CASE (3)
           
           IF ( this%restart_freq_time_wall <= 0.0_MK ) THEN
              PRINT *, "io_check_parameters : ", &
                   "restart_freq_time_wall <=0.0 !"
              io_check_parameters = .FALSE.
              GOTO 9999
           END IF
           
        END SELECT
        
9999    CONTINUE  
        
        RETURN        
        
      END FUNCTION io_check_parameters
