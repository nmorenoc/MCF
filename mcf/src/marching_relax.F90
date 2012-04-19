      SUBROUTINE marching_relax(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : marching_relax
        !----------------------------------------------------
        !
        ! Purpose     : Pre-run the simulation, in order
        !               to relax the particles' configuration.
        !
        ! Reference   : Initial ideas were stimulated by
        !               Sergej Litvinov.
        !
        ! Remark      : 1) Colloidal particles positions can
        !               be relaxed also, if the code has
        !               been compiled with a flag
        !               __COLLOID_RELAX
        !               2) If the fluid is non-Newtonian,
        !               we consider it Newtonian during relax
        !               run and switch it back to non-Newtonian
        !               afterwards.
        !
        !
        ! Revisions   : V0.1, 19.11.2010, small bug fixed for
        !               non-Newtonian relax run. Although it
        !               is considered Newtonian fluid during
        !               relax run, ct(conformation tensor)
        !               must be mapped to keep consistency.
        !
        !               V0.1 05.03 2010, original version.
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
        
        TYPE(Marching), INTENT(INOUT)   :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        
     	!----------------------------------------------------
      	! Local variables starts here :
      	!----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Pointer to Right hand side.
        ! Pointer to State Equation.
        !----------------------------------------------------
        
        TYPE(Rhs), POINTER              :: trhs
        TYPE(StateEquation), POINTER    :: tstateEquation
        
        !----------------------------------------------------
        ! Control parameters :
        !----------------------------------------------------
        
        LOGICAL                         :: symmetry_target
        LOGICAL                         :: symmetry
        LOGICAL                         :: dynamic_density_ref
        INTEGER                         :: stateEquation_type
        LOGICAL                         :: Newtonian
        LOGICAL                         :: Brownian
        LOGICAL                         :: p_energy
        INTEGER                         :: integrate_colloid_type
        
        !----------------------------------------------------
        ! Physics parameters :
        !----------------------------------------------------
        
        INTEGER                         :: num_species
        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(:), POINTER :: dx
        
        REAL(MK)                        :: dt
        
        REAL(MK)                        :: rho
        REAL(MK)			:: rho_min
        REAL(MK)			:: rho_max      
        REAL(MK)                        :: kt
        REAL(MK)                        :: c
        REAL(MK)                        :: rho_ref
        REAL(MK)                        :: gamma
        
        INTEGER                         :: relax_type
        REAL(MK)                        :: dt_relax
        INTEGER				:: step_current
        INTEGER                         :: step_relax
        REAL(MK)			:: time_current
        
        REAL(MK)                        :: disorder_level
        REAL(MK)                        :: disorder
        REAL(MK)                        :: kt_relax
        REAL(MK)                        :: c_relax
        
        REAL(MK)                        :: k_energy0
        REAL(MK)                        :: k_energy
        REAL(MK), DIMENSION(:), POINTER :: momentum
        
        !----------------------------------------------------
        ! body force type
        !----------------------------------------------------
        
        INTEGER                         :: body_force_type
        
     	!----------------------------------------------------
        ! Colloid parameters :
     	!----------------------------------------------------
        
        INTEGER                         :: num_colloid
        TYPE(Colloid), POINTER          :: colloids
        INTEGER                         :: coll_rho_type
        INTEGER                         :: coll_body_force_type
        LOGICAL                         :: coll_translate
        LOGICAL                         :: coll_rotate
        REAL(MK),DIMENSION(:,:,:),POINTER:: coll_v
        REAL(MK),DIMENSION(:,:,:),POINTER:: coll_omega
        REAL(MK), ALLOCATABLE, DIMENSION(:,:)   :: coll_drag
        REAL(MK), ALLOCATABLE, DIMENSION(:,:)   :: coll_torque
         REAL(MK),DIMENSION(:,:,:),POINTER:: coll_reset
        REAL(MK)                        :: coll_k
        REAL(MK), DIMENSION(3)          :: coll_mom        
        
        !----------------------------------------------------
        ! Boundary parameters :
     	!----------------------------------------------------
        
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: wall_rho_type
        INTEGER                         :: num_sym
        INTEGER                         :: num_wall_sym
        INTEGER                         :: num_wall_solid
        INTEGER                         :: num_le
        INTEGER                         :: num_shear
        REAL(MK),DIMENSION(:,:),POINTER :: shear_v0
        REAL(MK),DIMENSION(:,:),POINTER :: shear_reset;
        REAL(MK),DIMENSION(3,6)         :: wall_drag_c
        
        !----------------------------------------------------
        ! Number of real, all and ghost particles :
     	!----------------------------------------------------
        
        INTEGER                         :: num_part_real
        INTEGER                         :: num_part_all
        INTEGER                         :: num_part_ghost
        
        !----------------------------------------------------
        ! MPI parameters :
        !----------------------------------------------------
        
        INTEGER                         :: rank
        INTEGER                         :: comm
        INTEGER                         :: MPI_PREC
        
	!----------------------------------------------------
        ! writing frequency
     	!----------------------------------------------------

        INTEGER                         :: output_particles_freq_step
        INTEGER                         :: statistic_freq_step

        
        !----------------------------------------------------
        ! t_x : positions for particles, for building cell
        !       list.
        !----------------------------------------------------
        
        REAL(MK),DIMENSION(:,:),POINTER :: t_x
        

#ifdef __DEBUG
        
        !----------------------------------------------------
	! Debug variables.
        !----------------------------------------------------
        
        INTEGER                         :: debug_flag
        REAL(MK)			:: time_routine_start
        
#endif
        
        !----------------------------------------------------
        ! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info	= 0
        stat_info_sub	= 0
        
        NULLIFY(trhs)
        NULLIFY(tstateEquation)
        
        NULLIFY(dx)
        
        NULLIFY(momentum)
        
        num_colloid = 0
        NULLIFY(colloids)
        coll_rho_type = 0
        coll_translate = .FALSE.
        coll_rotate    = .FALSE.
        
        NULLIFY(coll_v)
        NULLIFY(coll_omega)
        NULLIFY(coll_reset)

        coll_k      = 0.0_MK
        coll_mom(:) = 0.0_MK
        
        NULLIFY(bcdef)
        NULLIFY(tboundary)
        wall_rho_type = 0
        NULLIFY(shear_v0)
        NULLIFY(shear_reset)
        wall_drag_c(:,:) = 0.0_MK
        
        NULLIFY(t_x)
        
        !----------------------------------------------------
        ! MPI parameters from the technique object.
        !----------------------------------------------------
        
        comm     = technique_get_comm(this%tech,stat_info_sub)
        rank     = technique_get_rank(this%tech,stat_info_sub)
        MPI_PREC = technique_get_MPI_PREC(this%tech,stat_info_sub)
        
        
#ifdef __DEBUG
        debug_flag = debug_get_flag(global_debug,stat_info_sub)
        IF( debug_flag == 3 ) THEN
           CALL debug_substart(global_debug,rank,&
                "marching_relax",&
                time_routine_start,stat_info_sub)
        END IF
#endif
        
        !----------------------------------------------------
        ! Control parameters :
        !
        ! Save target control parameters and set them 
        ! temporarily for relax run.
        !
        ! For relax run, we have Brownian + Newtonian fluid.
        ! Set symmetry to conserve momentum.
        !----------------------------------------------------
        
        symmetry_target     = &
             control_get_symmetry(this%ctrl,stat_info_sub)
        dynamic_density_ref = &
             control_get_dynamic_density_ref(this%ctrl,stat_info_sub)
        stateEquation_type   = &
             control_get_stateEquation_type(this%ctrl,stat_info_sub)
        Newtonian           = &
             control_get_Newtonian(this%ctrl,stat_info_sub)
        Brownian            = &
             control_get_Brownian(this%ctrl,stat_info_sub)
        p_energy            = &
             control_get_p_energy(this%ctrl,stat_info_sub)
        integrate_colloid_type = &
             control_get_integrate_colloid_type(this%ctrl,stat_info_sub)
      
        symmetry = .TRUE.
        CALL control_set_symmetry(this%ctrl,symmetry, stat_info_sub)
        CALL control_set_Brownian(this%ctrl,.TRUE., stat_info_sub)
        CALL control_set_Newtonian(this%ctrl,.TRUE., stat_info_sub)
        
        !----------------------------------------------------
        ! Get certain physics parameters.
        !
        ! Save desired dt, kt, c.
        !----------------------------------------------------
        
        num_species = &
             physics_get_num_species(this%phys,stat_info_sub)
        num_dim     = &
             physics_get_num_dim(this%phys,stat_info_sub)
        CALL physics_get_dx(this%phys,dx,stat_info_sub)
      
        dt      = &
             physics_get_dt(this%phys,stat_info_sub)
        rho     = &
             physics_get_rho(this%phys,stat_info_sub)
        kt      = &
             physics_get_kt(this%phys,stat_info_sub)
        c       = &
             physics_get_c(this%phys,stat_info_sub)
        rho_ref = &
             physics_get_rho_ref(this%phys,stat_info_sub)
        gamma   =  &
             physics_get_gamma(this%phys,stat_info_sub)

        !----------------------------------------------------
        ! Relax parameters :
        !----------------------------------------------------
        
        relax_type     = &
             physics_get_relax_type(this%phys,stat_info_sub)
        dt_relax       = &
             physics_get_dt_relax(this%phys,stat_info_sub)
        step_current   = 0
        step_relax     = &
             physics_get_step_relax(this%phys,stat_info_sub)
        time_current   = 0.0_MK
        disorder_level = &
             physics_get_disorder_level(this%phys,stat_info_sub)
       
        kt_relax       = &
             physics_get_kt_relax(this%phys,stat_info_sub)
      
        c_relax        = &
             physics_get_c_relax(this%phys,stat_info_sub)
        
        !----------------------------------------------------
        ! 1) Set state Equation according to relax run
        ! sound speed.
        !----------------------------------------------------
        
        CALL particles_get_stateEquation(this%particles,&
             tstateEquation,stat_info_sub)
        CALL stateEquation_new(tstateEquation,&
             stateEquation_type,&
             c_relax,rho,rho_ref,gamma,stat_info_sub)
        
        !----------------------------------------------------
        ! 2) Set dt and kt in right hand side obeject.
        ! where dt is used for random force.
        ! Since dt can be a variable during simulation,
        ! (e.g, dt for relax run is different for real run,
        !  or adaptive time step)
        ! it is necessary to update dt in rhs.
        !
        ! 3)Set kt, since kt can be different for relax
        ! run simulation and real run simulation.
        !----------------------------------------------------
         
        CALL particles_get_rhs(this%particles,trhs,stat_info_sub)
        CALL rhs_set_Brownian(trhs,.TRUE.,stat_info_sub)
        CALL rhs_set_dt(trhs,dt_relax,stat_info_sub)
        CALL rhs_set_kt(trhs,kt_relax,stat_info_sub)
        
        !----------------------------------------------------
        ! Body force parameters :
        ! Set body force type zero, i.e., no body force.
        !----------------------------------------------------
     
        body_force_type = &
             physics_get_body_force_type(this%phys,stat_info_sub)
        CALL physics_set_body_force_type(this%phys,0,stat_info_sub)
        
        !----------------------------------------------------
        ! Colloid parameters :
        !
        ! For complex fluids, return a pointer of colloid.
        !----------------------------------------------------
        
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF ( num_species > 1 .AND. num_colloid > 0 ) THEN
           
           CALL physics_get_colloid(this%phys,&
                colloids,stat_info_sub)
           
           coll_rho_type = &
                colloid_get_rho_type(colloids,stat_info_sub)
           
           coll_body_force_type = &
                colloid_get_body_force_type(colloids,stat_info_sub)
         
           !-------------------------------------------------
           ! Save desirded status of translate and rotate
           ! of colloids, reset them to non-moveable,
           ! prepare for relax fluid particles.
           !-------------------------------------------------
           
           coll_translate = &
                colloid_get_translate(colloids,stat_info_sub)
           coll_rotate    = &
                colloid_get_rotate(colloids,stat_info_sub)
           CALL colloid_get_v(colloids, &
                coll_v,stat_info_sub)
           CALL colloid_get_omega(colloids, &
                coll_omega,stat_info_sub)
           
           CALL colloid_set_translate(colloids, &
                .FALSE.,stat_info_sub)
           CALL colloid_set_rotate(colloids, &
                .FALSE.,stat_info_sub)
           
           CALL colloid_set_body_force_type(colloids, &
                0,stat_info_sub)
          
#if __COLLOID_RELAX
           CALL colloid_set_translate(colloids, &
                .TRUE.,stat_info_sub)
           CALL colloid_set_rotate(colloids, &
                .TRUE.,stat_info_sub)
#endif
           
           ALLOCATE(coll_reset(3,num_colloid,integrate_colloid_type))
           coll_reset(:,:,:) = 0.0_MK
           CALL colloid_set_v(colloids, &
                coll_reset(1:num_dim,1:num_colloid,1:integrate_colloid_type),&
                stat_info_sub)
           CALL colloid_set_omega(colloids, &
                coll_reset(1:3,1:num_colloid,1:integrate_colloid_type),&
                stat_info_sub)
           
           ALLOCATE(coll_drag(num_dim,num_colloid))
           ALLOCATE(coll_torque(3,num_colloid))
      
        END IF ! num_colloid > 0
        
        
        !----------------------------------------------------
        ! Boundary parameters :
        !----------------------------------------------------
        
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        
        !----------------------------------------------------
        ! Save desired shear_v0 and reset it to zero for
        ! relaxing particles' configuration.
        !----------------------------------------------------
        
        CALL boundary_get_shear_v0(tboundary, &
             shear_v0,stat_info_sub)
        ALLOCATE(shear_reset(num_dim,2*num_dim))
        shear_reset(:,:) = 0.0_MK
        CALL boundary_set_shear_v0(tboundary, &
             shear_reset,stat_info_sub)
        
        num_sym      = &
             boundary_get_num_sym(tboundary,stat_info_sub)
        num_wall_sym = &
             boundary_get_num_wall_sym(tboundary,stat_info_sub)
        num_wall_solid = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        
        IF (num_wall_solid > 0) THEN
           
           wall_rho_type =  &
                boundary_get_wall_rho_type(tboundary,stat_info_sub)
           
        END IF
        
        num_le       = &
             boundary_get_num_le(tboundary,stat_info_sub)
        num_shear    = &
             boundary_get_num_shear(tboundary,stat_info_sub)
        
        !----------------------------------------------------
        ! Update boundary :
        !
        ! Walls' velocity in case of oscillating shear;
        ! Sheared length, in case of Lees-Edwards boundary.
        !
        ! Reset them all to values at time zero.
        !----------------------------------------------------
        
        IF ( num_shear > 0 ) THEN
           
           CALL boundary_update_boundary(tboundary,&
                0.0_MK,0.0_MK,stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_relax: ", &
                   "boundary updating boundary failed!"
              stat_info = -1           
              GOTO 9999           
           END IF
           
        END IF
        
        output_particles_freq_step = &
             io_get_output_particles_relax_freq_step(this%io,stat_info_sub)
        statistic_freq_step = &
             io_get_statistic_relax_freq_step(this%io,stat_info_sub)
  
        !----------------------------------------------------
      	! Create ghost particles on righ-inner
        ! boundaries (symmetry) for each process.
        ! Ghosts are with position, mass and IDs.
        !
        ! Even using single process, this has to be done,
        ! since it guarantees the boundary condition.
        !
	! Has to be done BEFORE building neighbor list.
      	!----------------------------------------------------
        
        CALL particles_map_ghost_get(this%particles, &
             l_map_x  = .TRUE., l_map_m = .TRUE., &
             l_map_id = .TRUE., stat_info=stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "marching_relax: ", &
                "Creating ghosts with x, m, IDs failed!"
           stat_info = -1           
           GOTO 9999           
        END IF
        
        !----------------------------------------------------
        ! According to different boundary condition,
        ! we have set up the IDs of boundary particles
        ! which are ghosts also.
        !----------------------------------------------------
        
        CALL particles_set_boundary_ghost_id(this%particles, &
             stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "marching_relax: ", &
                "Setting boundary ghosts ID failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! The number of real particles on this process.
        ! It should not be changed yet;
        ! The number of ghost particles on this process has
        ! been calculated newly.
        ! The number of all particles of this process has 
        ! been calculated newly also, however,
        ! num_part_all = num_part_real + num_part_ghost
        ! is always true.
        !----------------------------------------------------
        
        num_part_real  = &
             particles_get_num_part_real(this%particles,stat_info_sub)
        num_part_ghost = &
             particles_get_num_part_ghost(this%particles,stat_info_sub)
        num_part_all  = &
             particles_get_num_part_all(this%particles,stat_info_sub)
     
        
        !----------------------------------------------------
        ! Get all particles' positions (including ghosts),
      	! to build neighbor list.
      	!----------------------------------------------------
        
        CALL particles_get_x(this%particles,t_x, &
             num_part_all,stat_info_sub)        
        CALL technique_build_list(this%tech,t_x, &
             num_part_all,symmetry,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *,"marching_relax: ", &
                "Building neighbour list failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! Compute mass density/number density at time_current.
        !----------------------------------------------------
        
        CALL particles_compute_density(this%particles, &
             stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN           
           PRINT *, "marching_relax: ", &
                "Computing density failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! For symmtery inter-process communication :
        !
        ! swap the send and receive buffer,in order to send
        ! contribution of ghost particles to their hosts
        ! on neigboring processes and receive contribution 
        ! of ghosts from other neigboring processes.
        !
        ! Some local particles which don't have contribution
        ! from other processes, will stay the same.
        !
        ! Even using single process, this has to be done,
        ! since it guarantees the boundary conditions also.
      	!----------------------------------------------------
        
        IF ( symmetry ) THEN
           
           CALL particles_map_ghost_put(this%particles,&
                l_map_x = .TRUE., l_map_rho = .TRUE., &
                stat_info=stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN           
              PRINT *, "marching_relax: ",&
                   "Receiveing contribution from ghosts failed!"
              stat_info = -1
              GOTO 9999           
           END IF
           
        END IF
        
        !-------------------------------------------------
        ! Set colloidal boundary particles velocity
        ! according to translation and rotation speed.
        !-------------------------------------------------
        
        IF ( num_colloid > 0 ) THEN
           
           CALL particles_set_colloid_velocity(this%particles,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'marching_relax: ',&
                   'Setting colloid velocity failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        IF ( num_shear > 0 ) THEN
           
           CALL particles_set_boundary_velocity(this%particles,&
                stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_relax: ", &
                   "particles setting boundary failed!"
              stat_info = -1           
              GOTO 9999           
           END IF

        END IF
        
        !----------------------------------------------------
        ! After the densities are updated on each proceses, 
        ! the values have to be sent to ghost particles at 
        ! neigboring processes.
        !
        ! At this point, also update ghosts with velocity,
        ! which is needed for vgt(Oldroyd-B model) and 
        ! interaction(force etc.) calculation.
        !
        ! And update ghosts with ct for calculating force in 
        ! case of non-Newtonian viscoelastic Oldroyd-B model.
        !
        ! Even using single process, this has to be done,
        ! since it guarantee the boundary conditions also.
     	!----------------------------------------------------
        
        CALL particles_map_ghost_get(this%particles,&
             l_map_x  = .TRUE., l_map_rho=.TRUE., &
             l_map_v  = .TRUE., &
             l_map_ct = (.NOT. Newtonian), &
             stat_info=stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN           
           PRINT *, 'marching_relax: ', &
                'Updating ghosts with x, rho, v failed!'
           stat_info = -1           
           GOTO 9999           
        END IF
        
        
        !----------------------------------------------------
        ! The number of real, all and ghost particles
        ! should not be changed, shown here for clarity.
        ! num_part_all = num_part_real + num_part_ghost
        !----------------------------------------------------
        
        num_part_real  = &
             particles_get_num_part_real(this%particles,stat_info_sub)
        num_part_ghost = &
             particles_get_num_part_ghost(this%particles,stat_info_sub)
        num_part_all  = &
             particles_get_num_part_all(this%particles,stat_info_sub)
        
        !----------------------------------------------------
        ! After mapping velocity, according to different 
        ! boundary conditions, we have to set up velocity 
        ! of boundary particles which are ghosts also.
        !----------------------------------------------------
        
        CALL particles_set_boundary_ghost_velocity(this%particles, &
             stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN           
           PRINT *, "marching_relax : ", &
                "Setting boundary ghosts velocity failed!"
           stat_info = -1
           GOTO 9999           
        END IF
        
        !----------------------------------------------------
        ! If reference density rho_ref given from input file is
        ! negative, we find the minimum density during
        ! simulation and set by particles_set_rho_ref
        ! into state equation.
        !----------------------------------------------------
        
        IF ( dynamic_density_ref ) THEN
           
           CALL particles_find_density_extreme(this%particles, &
                comm, MPI_PREC, stat_info_sub)
           
           IF ( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ", &
                   "Finding density extrem failed!"
              stat_info = -1
              GOTO 9999
           END IF

           CALL particles_set_stateEquation_rho_ref(this%particles, &
                stat_info_sub)
           
           IF ( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ", &
                   "Setting density extreme failed!"
              stat_info = -1
              GOTO 9999
           END IF

        END IF


        !----------------------------------------------------
      	! Compute pressure (including ghost particles)
        ! at time_current, since the pressure of ghost particles 
        ! are used for calculating the interaction(force) also.
        !----------------------------------------------------
        
        CALL particles_compute_pressure(this%particles,&
             num_part_all,stat_info_sub)
        
        IF ( stat_info_sub /=0 ) THEN
           PRINT *, "marching_relax: ", &
                "Computing pressure failed!"
           stat_info = -1
           GOTO 9999     
        END IF
        
        !----------------------------------------------------
      	! Compute interactions between particles 
        ! at time_current.
        ! e.g., force.
      	!----------------------------------------------------
        
        CALL particles_compute_interaction(this%particles, &
             stat_info_sub)
        
        IF ( stat_info_sub /=0 ) THEN
           PRINT *, "marching_relax: ", &
                "Computing interaction failed!"
           stat_info = -1
           GOTO 9999
        END IF

        
        !----------------------------------------------------
        ! For symmtery inter-process communication:
        !----------------------------------------------------
        
        IF( symmetry ) THEN
           
           !-------------------------------------------------
           ! Swap the send and receive buffer, in order to 
           ! send contribution of force on ghost particles 
           ! to their host processes and receive contribution
           ! of force from other processes.
           !
           ! Pontential energy depends on requirement.
           !
           ! Some particles which don't have contribution 
           ! from other processes, remain the same.
           !
           ! Even using single process, this has to be done,
           ! since it guarantee the boundary condition also.
           !-------------------------------------------------
           
!#ifdef __PARTICLES_FORCE_SEPARATE
!           CALL particles_map_ghost_put(this%particles, &
!                l_map_x   = .TRUE., l_map_f  = .TRUE., &
!                l_map_fp  = .TRUE., l_map_fv = .TRUE., &
!                l_map_au  = p_energy, &
!                stat_info = stat_info_sub)           
!#else
           CALL particles_map_ghost_put(this%particles, &
                l_map_x   = .TRUE., l_map_f  = .TRUE., &
                l_map_au  = p_energy, &
                stat_info = stat_info_sub)
!#endif
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, 'marching_relax: ',&
                   'Receiving force(au) from ghosts failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
        !----------------------------------------------------
        ! Since paire-wise force have been calculated,
        ! sum up the force on colloids,
        ! if colloids are present.
        ! Compute colloid-colloid and colloid-wall 
        ! interactions.
        !----------------------------------------------------
     
        IF ( num_colloid > 0 ) THEN
           
           !-------------------------------------------------
           ! Sum up force/torque exerted on parts of
           ! colloids on this process.
           !-------------------------------------------------
           
           CALL particles_collect_colloid_interaction(this%particles,&
                coll_drag,coll_torque,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ",&
                   "Summing up interaction on colloid locally has problem!"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Sum up force/torque exerted on colloids
           ! from all processes.
           !-------------------------------------------------
           
           CALL colloid_collect_particles_interaction(colloids,&
                comm,MPI_PREC,coll_drag,coll_torque,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ",&
                   "Summing up interaction on colloid globally has problem!"
              stat_info = -1
              GOTO 9999
           END IF
   
           !-------------------------------------------------
           ! Compute colloid-colloid and colloid-wall
           ! interactions.
           !-------------------------------------------------
           
           CALL colloid_compute_interaction(colloids,comm,&
                MPI_PREC, coll_drag, coll_torque,&
                wall_drag_c(1:num_dim,1:num_dim*2),&
                stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ",&
                   "c-c or c-w interaction has problem!"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Compute colloids accelerations, i.e.,
           ! translation and rotation.
           !-------------------------------------------------
           
           CALL colloid_compute_acceleration(colloids,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ",&
                   "Computing colloids accelerations has problem!"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! num_colloid > 0

        !----------------------------------------------------
        ! If there is wall using symmetry or solid boundary
        ! particles, or Lees-Edwards boundary particles.
        ! we reset boundary interaction to zero.
        !----------------------------------------------------
        
        IF ( num_shear > 0 ) THEN
           
           CALL particles_reset_boundary_interaction(this%particles, &
                stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ",&
                   "Resetting boundary particles interaction failed!"
              stat_info = -1
              GOTO 9999
           END IF

           CALL particles_reset_boundary_ghost_interaction(this%particles, &
                stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ",&
                   "Resetting boundary ghost particles interaction failed!"
              stat_info = -1
              GOTO 9999
           END IF
        
        END IF
        
        
        CALL statistic_compute_disorder(this%statis, &
             this%particles,dx,stat_info_sub)
        
        disorder =  &
             statistic_get_disorder_square_root(this%statis, &
             stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN              
           PRINT *, "marching_relax: ",&
                "Computing disorder failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        CALL statistic_compute_statistic(this%statis, &
             this%particles,coll_k,coll_mom(1:num_dim), &
             stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN              
           PRINT *, "marching_relax: ",&
                "Computing statistic failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( rank == 0 ) THEN
           
           CALL io_open_statistic_relax(this%io,rank,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ",&
                   "Opening statistic_relax file failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
           CALL io_write_statistic_relax(this%io,rank, &
                step_current,time_current, &
                this%statis,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ",&
                   "Writting statistic_relax file failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        num_part_real = &
             particles_get_num_part_real(this%particles,stat_info_sub)
        

        CALL io_write_particles_relax(this%io,rank, &
             step_current,this%particles, &
             num_part_real, stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN
           PRINT *, "marching_relax: ",&
                "Writting particles_relax file failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        !----------------------------------------------------
        ! ###################################################
        ! ### Time loop for integration for relaxation.   ###
        ! ###################################################
        !----------------------------------------------------
      	!----------------------------------------------------
	
        DO WHILE ( ( relax_type == 1 .AND. &
             step_current <  step_relax )  .OR.  &
             ( relax_type == 2 .AND. disorder > disorder_level ) )
           
           CALL marching_integrate(this,step_current,&
                time_current,dt_relax,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN              
              PRINT *, "marching_relax: ",&
                   "Integrating failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Update current step and time.
           !-------------------------------------------------

           step_current = step_current + 1
           time_current = time_current + dt_relax
           
           IF ( MOD(step_current,output_particles_freq_step) == 0 ) THEN
              
              num_part_real = &
                   particles_get_num_part_real(this%particles,stat_info_sub)
              
              CALL io_write_particles_relax(this%io,rank, &
                   step_current,this%particles, &
                   num_part_real, stat_info_sub)
              
              IF ( stat_info_sub /=0 ) THEN
                 PRINT *, "marching_relax: ",&
                      "Writting particles_relax file failed!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF
           
           CALL statistic_compute_disorder(this%statis, &
                this%particles,dx,stat_info_sub)
           
           disorder =  &
                statistic_get_disorder_square_root(this%statis, &
                stat_info_sub)
           
           IF ( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ",&
                   "Computing disorder failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
           IF ( MOD(step_current,statistic_freq_step) == 0 ) THEN
              
              CALL statistic_compute_statistic(this%statis, &
                   this%particles,coll_k,coll_mom(1:num_dim), &
                   stat_info_sub)
              
              IF ( stat_info_sub /=0 ) THEN              
                 PRINT *, "marching_relax: ",&
                      "Computing statistic failed!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              IF ( dynamic_density_ref ) THEN
                 
                 rho_min = particles_get_rho_min(this%particles, &
                      stat_info_sub)
                 rho_max = particles_get_rho_max(this%particles, &
                      stat_info_sub)
                 
                 CALL statistic_set_rho_min(this%statis,rho_min, &
                      stat_info_sub)
                 
                 CALL statistic_set_rho_max(this%statis,rho_max, &
                      stat_info_sub)
           
              END IF
              
              IF ( rank == 0 ) THEN
                 
                 CALL io_write_statistic_relax(this%io,rank, &
                      step_current,time_current, &
                      this%statis,stat_info_sub)
                 
                 IF( stat_info_sub /=0 ) THEN
                    PRINT *, "marching_relax : ",&
                         "Writting statistic_relax file failed !"
                    stat_info = -1
                    GOTO 9999
                 END IF
              
              END IF
              
           END IF
           
        END DO ! step_current < step_relax
        
        !----------------------------------------------------
        !----------------------------------------------------
        !####################################################
        !#######   End of time loop for integration.  #######
        !####################################################
        !----------------------------------------------------
        !----------------------------------------------------
        
        CALL statistic_get_statistic(this%statis,&
             k_energy,momentum,stat_info_sub)
        
        k_energy0 = k_energy
        
        CALL control_set_Brownian(this%ctrl,.FALSE., stat_info_sub)
        CALL rhs_set_Brownian(trhs,.FALSE.,stat_info_sub)
        !CALL rhs_set_dt(trhs,dt,stat_info_sub)
        !CALL rhs_set_kt(trhs,0.0_MK,stat_info_sub)
      
        !----------------------------------------------------
        !----------------------------------------------------
        ! ###################################################
        ! ### Integration for reducing kinetic energy.    ###
        ! ###################################################
        !----------------------------------------------------
      	!----------------------------------------------------
	
        DO WHILE ( k_energy * mcf_kinetics_reduce_factor > k_energy0  )
           
           CALL marching_integrate(this,step_current,&
                time_current,dt_relax,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN              
              PRINT *, "marching_relax : ",&
                   "Integrating failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Update current step and time.
           !-------------------------------------------------

           step_current = step_current + 1
           time_current = time_current + dt_relax
           
           IF ( MOD(step_current,output_particles_freq_step) == 0 ) THEN
              
              num_part_real = &
                   particles_get_num_part_real(this%particles,stat_info_sub)
              
              CALL io_write_particles_relax(this%io,rank, &
                   step_current,this%particles, &
                   num_part_real, stat_info_sub)
              
              IF ( stat_info_sub /=0 ) THEN
                 PRINT *, "marching_relax: ",&
                      "Writting particles_relax file failed!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF
           
           CALL statistic_compute_disorder(this%statis, &
                this%particles,dx,stat_info_sub)
           
           disorder =  &
                statistic_get_disorder_square_root(this%statis, &
                stat_info_sub)
           
           IF ( stat_info_sub /=0 ) THEN              
              PRINT *, "marching_relax: ",&
                   "Computing disorder failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
           IF ( MOD(step_current,statistic_freq_step) == 0 ) THEN
              
              CALL statistic_compute_statistic(this%statis, &
                   this%particles,coll_k,coll_mom(1:num_dim), &
                   stat_info_sub)
              
              IF ( stat_info_sub /=0 ) THEN              
                 PRINT *, "marching_relax: ",&
                      "Computing statistic failed!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              CALL statistic_get_statistic(this%statis,&
                   k_energy,momentum,stat_info_sub)
              
              IF ( dynamic_density_ref ) THEN
                 
                 rho_min = particles_get_rho_min(this%particles, &
                      stat_info_sub)
                 rho_max = particles_get_rho_max(this%particles, &
                      stat_info_sub)
                 
                 CALL statistic_set_rho_min(this%statis,rho_min, &
                      stat_info_sub)
                 
                 CALL statistic_set_rho_max(this%statis,rho_max, &
                      stat_info_sub)
           
              END IF
              
              IF ( rank == 0 ) THEN
                 
                 CALL io_write_statistic_relax(this%io,rank, &
                      step_current,time_current, &
                      this%statis,stat_info_sub)
                 
                 IF ( stat_info_sub /=0 ) THEN
                    PRINT *, "marching_relax: ",&
                         "Writting statistic_relax file failed!"
                    stat_info = -1
                    GOTO 9999
                 END IF
              
              END IF
              
           END IF
           
        END DO ! step_current < step_relax
        
        !----------------------------------------------------
        !----------------------------------------------------
        !####################################################
        !#######   End of time loop for integration.  #######
        !####################################################
        !----------------------------------------------------
        !----------------------------------------------------
        
        IF ( MOD(step_current,output_particles_freq_step) /= 0 ) THEN
           
           num_part_real = &
                particles_get_num_part_real(this%particles,stat_info_sub)
           
           CALL io_write_particles_relax(this%io,rank, &
                step_current,this%particles, &
                num_part_real, stat_info_sub)
           
           IF ( stat_info_sub /=0 ) THEN
              PRINT *, "marching_relax: ",&
                   "Writting particles_relax file failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
        IF ( MOD(step_current,statistic_freq_step) /= 0 ) THEN
           
           CALL statistic_compute_statistic(this%statis, &
                this%particles, coll_k,coll_mom(1:num_dim), &
                stat_info_sub)
           
           IF ( stat_info_sub /=0 ) THEN              
              PRINT *, "marching_relax: ",&
                   "Computing statistic failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
           CALL statistic_compute_disorder(this%statis, &
                this%particles,dx,stat_info_sub)
           
           IF ( stat_info_sub /=0 ) THEN              
              PRINT *, "marching_relax: ",&
                   "Computing disorder failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
           IF ( dynamic_density_ref ) THEN
              
              rho_min = particles_get_rho_min(this%particles, &
                   stat_info_sub)
              rho_max = particles_get_rho_max(this%particles, &
                   stat_info_sub)
              
              CALL statistic_set_rho_min(this%statis,rho_min, &
                   stat_info_sub)
              
              CALL statistic_set_rho_max(this%statis,rho_max, &
                   stat_info_sub)
              
           END IF
           
           IF ( rank == 0 ) THEN
              
              CALL io_write_statistic_relax(this%io,rank, &
                   step_current,time_current, &
                   this%statis,stat_info_sub)
              
              IF ( stat_info_sub /=0 ) THEN
                 PRINT *, "marching_relax: ",&
                      "Writting statistic_relax file failed!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF
           
        END IF


        !----------------------------------------------------
        ! Reset particles velocity and force to zero
        !----------------------------------------------------
        
        CALL particles_reset_v(this%particles,stat_info_sub)
        CALL particles_reset_f(this%particles,stat_info_sub)
        
        !----------------------------------------------------
        ! Write restart file after relax.
        !----------------------------------------------------
        
        num_part_real = &
             particles_get_num_part_real(this%particles,stat_info_sub)
        
        CALL io_write_restart_particles_relax(this%io,rank, &
             step_current,this%particles, &
             num_part_real, stat_info_sub)
        
        IF ( stat_info_sub /=0 ) THEN
           PRINT *, "marching_relax: ",&
                "Writting restart_particles_relax file failed!"
           stat_info = -1
           GOTO 9999
        END IF
           
        !----------------------------------------------------
        ! Restore parameters to prepare for real running.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Restore control parameters.
        !----------------------------------------------------

        CALL control_set_symmetry(this%ctrl,symmetry_target,stat_info_sub)
        CALL control_set_Brownian(this%ctrl,Brownian,stat_info_sub)
        CALL control_set_Newtonian(this%ctrl,Newtonian,stat_info_sub)
        
        !----------------------------------------------------
        ! Restore state equation using real sound speed.
        !----------------------------------------------------

        CALL stateEquation_new(tstateEquation,&
             stateEquation_type,&
             c,rho,rho_ref,gamma,stat_info_sub)
        
        !----------------------------------------------------
        ! Restore RHS variables.
        !----------------------------------------------------
        
        CALL rhs_set_Brownian(trhs,Brownian,stat_info_sub)
        CALL rhs_set_dt(trhs,dt,stat_info_sub)
        CALL rhs_set_kt(trhs,kt,stat_info_sub)
       
        !----------------------------------------------------
        ! Restore body force type.
        !----------------------------------------------------
   
        CALL physics_set_body_force_type(this%phys, &
             body_force_type,stat_info_sub)
        
        !----------------------------------------------------
        ! Restore shear velocity.
        !----------------------------------------------------
   
        CALL boundary_set_shear_v0(tboundary, &
             shear_v0, stat_info_sub)
        
        !----------------------------------------------------
        ! Restore colloid properties.
        !----------------------------------------------------

        IF (  num_species > 1 .AND. num_colloid > 0 ) THEN
           
           CALL colloid_set_translate(colloids, &
                coll_translate,stat_info_sub)
           CALL colloid_set_rotate(colloids, &
                coll_rotate,stat_info_sub)
           CALL colloid_set_body_force_type(colloids, &
                coll_body_force_type,stat_info_sub)        
           CALL colloid_set_v(colloids, &
                coll_v,stat_info_sub)
           CALL colloid_set_omega(colloids, &
                coll_omega,stat_info_sub)
           
        END IF ! num_species > 1
        
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Close statistic file.
        !----------------------------------------------------
        
        CALL io_close_statistic_relax(this%io,stat_info_sub)
        
        IF (stat_info_sub /=0 ) THEN
           PRINT *,"marching_marching: ",&
                "Closing statistic_relax files failed!"
           stat_info = -1
        END IF
        
        !----------------------------------------------------
        ! Release all dynamics memories.
        !----------------------------------------------------
        
        IF(ASSOCIATED(dx)) THEN
           DEALLOCATE(dx)
        END IF
        
        IF(ASSOCIATED(momentum)) THEN
           DEALLOCATE(momentum)
        END IF
        
        IF(ASSOCIATED(coll_v)) THEN
           DEALLOCATE(coll_v)
        END IF

        IF(ASSOCIATED(coll_omega)) THEN
           DEALLOCATE(coll_omega)
        END IF
        
        IF(ASSOCIATED(coll_reset)) THEN
           DEALLOCATE(coll_reset)
        END IF
        
        IF(ASSOCIATED(t_x)) THEN
           DEALLOCATE(t_x)
        END IF
        
        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF
        
        IF(ASSOCIATED(shear_v0)) THEN
           DEALLOCATE(shear_v0)
        END IF
        
        IF(ASSOCIATED(shear_reset)) THEN
           DEALLOCATE(shear_reset)
        END IF
   
#ifdef __DEBUG        
        IF(debug_flag == 3 ) THEN
           CALL debug_substop(global_debug,rank, &
                'marching_relax:', &
                time_routine_start, stat_info_sub)
        END IF 
#endif
        
        RETURN
        
      END SUBROUTINE marching_relax
      
      
      
