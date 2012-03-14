      SUBROUTINE marching_marching(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : marching_marching
        !----------------------------------------------------
        !
        ! Purpose     : Time-Integration preparing, controling
        !               and finalizing routine.
        !
        !               Initial calculation of density, 
        !               interaction(force, etc) and iteration/
        !               evolution of step/time for integrating
        !               in time are done here; 
        !               Initial output and final output are 
        !               done here as well.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.5 07.12 2009, merge general
        !               boundary conditions, i.e., periodic,
        !               symmetry, wall using symmetry, solid
        !               wall and Lees-Edwards boundaries.
        !
        !               V0.4 30.07 2009
        !               Add computing velocity gradient 
        !               tensor which is used for Oldroyd-B 
        !               model of Viscoelastic Non-Newtonian 
        !               fluid.
        !
        !               V0.3 23.07 2009,
        !               Merge marching_marching_nonsym() and
        !               marching_marching_sym() together
        !
        !               V0.2 09.07 2009, 
        !               check again the work flow is correct
        !               and supply with more comments.
        !
        !               V0.1 15.07 2009, original version.
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
        ! Control parameters :
        !
     	!----------------------------------------------------
        
        LOGICAL                         :: read_external
        LOGICAL                         :: dynamic_density_ref
        LOGICAL                         :: symmetry
        LOGICAL                         :: Newtonian
        LOGICAL                         :: Brownian
        LOGICAL                         :: stress_tensor
        LOGICAL                         :: stress_tensor_p
        LOGICAL                         :: stress_tensor_v
        LOGICAL                         :: stress_tensor_r
        LOGICAL                         :: p_energy
        LOGICAL                         :: flow_v_fixed
        INTEGER                         :: adaptive_dt

        !----------------------------------------------------
        ! Physics parameters :
        !
    	!----------------------------------------------------
        
        INTEGER                         :: num_species
        INTEGER                         :: num_dim
        REAL(MK)                        :: dt
        INTEGER                         :: step_start
        INTEGER				:: step_current
        INTEGER                         :: step_end
        REAL(MK)			:: time_current
        REAL(MK)			:: time_end
        REAL(MK)			:: rho_min
        REAL(MK)			:: rho_max
        REAL(MK)			:: dt_f, dt_f_c
        
        !----------------------------------------------------
        ! non-Newtonian viscoelastic Oldroyd-B model
        ! parameters :
        !----------------------------------------------------
        
        LOGICAL                         :: eigen_dynamics
        
        !----------------------------------------------------
        ! Right hand side pointer
        !----------------------------------------------------

        TYPE(Rhs), POINTER              :: trhs
        
        !----------------------------------------------------
        ! body force and flow condition.
        !----------------------------------------------------
        
        REAL(MK), DIMENSION(:), POINTER :: body_force
        REAL(MK), DIMENSION(:), POINTER :: body_force_d
        INTEGER                         :: flow_adjust_freq
        
        !----------------------------------------------------
        ! Colloid parameters :
        !
        ! coll_k   : total kinetic energy of colloids
        ! coll_mom : total momentum of collods.
     	!----------------------------------------------------
        
        INTEGER                         :: num_colloid
        TYPE(Colloid), POINTER          :: colloids
        REAL(MK),DIMENSION(:,:),POINTER :: coll_drag
        REAL(MK),DIMENSION(:,:),POINTER :: coll_torque
        
        REAL(MK)                        :: coll_k
        REAL(MK), DIMENSION(:), POINTER :: coll_mom        
        
        !----------------------------------------------------
        ! Boundary parameters :
        ! wall_drag_p : drag from SDPD/SPH particles
        !               (solvent, colloidal boundary particles)
        ! wall_drag_c : drag from colloidal particles.
        !
        ! wall_drag_pp: pressure force from SDPD/SPH particles.
        ! wall_drag_pv: viscous force from SDPD/SPH particles.
        ! wall_drag_pr: random force from SDPD/SPH particles.
     	!----------------------------------------------------
        
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: num_sym
        INTEGER                         :: num_wall_sym
        INTEGER                         :: num_wall_solid
        INTEGER                         :: num_le
        INTEGER                         :: num_shear
        REAL(MK),DIMENSION(3,6)         :: wall_drag_p
        REAL(MK),DIMENSION(3,6)         :: wall_drag_c
#ifdef __WALL_FORCE_SEPARATE
        REAL(MK),DIMENSION(3,6)         :: wall_drag_pp
        REAL(MK),DIMENSION(3,6)         :: wall_drag_pv
        REAL(MK),DIMENSION(3,6)         :: wall_drag_pr
#endif

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
        ! Output condition flags :
     	!----------------------------------------------------
        
        LOGICAL                         :: l_write_particles
        LOGICAL                         :: l_write_conformation
        LOGICAL                         :: l_write_colloid
        LOGICAL                         :: l_write_statistic
        LOGICAL                         :: l_write_boundary
        LOGICAL                         :: l_write_restart_physics
        LOGICAL                         :: l_write_restart_particles
        LOGICAL                         :: l_write_restart_conformation
        
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
      
        REAL(MK),DIMENSION(:),POINTER   :: t_rho
        REAL(MK),DIMENSION(:,:),POINTER :: t_f
        INTEGER,DIMENSION(:,:),POINTER  :: t_id
#endif
        
        !----------------------------------------------------
        ! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info	= 0
        stat_info_sub	= 0
        
        NULLIFY(trhs)
        
        NULLIFY(body_force)
        NULLIFY(body_force_d)
       
        num_colloid = 0
        NULLIFY(colloids)
        NULLIFY(coll_drag)
        NULLIFY(coll_torque)
        coll_k  = 0.0_MK
        NULLIFY(coll_mom)
        
        NULLIFY(bcdef)
        NULLIFY(tboundary)
        wall_drag_p(:,:) = 0.0_MK
        wall_drag_c(:,:) = 0.0_MK
#ifdef __WALL_FORCE_SEPARATE
        wall_drag_pp(:,:) = 0.0_MK
        wall_drag_pv(:,:) = 0.0_MK
        wall_drag_pr(:,:) = 0.0_MK
#endif
        
        NULLIFY(t_x)
        
        !----------------------------------------------------
        ! MPI parameters from the technique object.
        !----------------------------------------------------
        
        comm     = technique_get_comm(this%tech,stat_info_sub)
        rank     = technique_get_rank(this%tech,stat_info_sub)
        MPI_PREC = technique_get_MPI_PREC(this%tech,stat_info_sub)
        
        
#ifdef __DEBUG
        debug_flag = debug_get_flag(global_debug,stat_info_sub)
        IF(debug_flag == 2 ) THEN
           CALL debug_substart(global_debug,rank,&
                "marching_marching",&
                time_routine_start,stat_info_sub)
        END IF

        NULLIFY(t_rho)
        NULLIFY(t_f)
        NULLIFY(t_id)
#endif
        
        !----------------------------------------------------
        ! Control parameters :
        !
        ! read_external   : reading particles from external
        !                   files or not.
        ! symmetry        : symmetry or not for inter-process
        !                   communication.
        ! Newtonian       : Newtonian or non-Newtonian fluid.
        ! p_energy        : wether we need potential energy.
        ! flow_v_fixed    : fixed flow velocity or not.
        !----------------------------------------------------
        
        read_external = &
             control_get_read_external(this%ctrl,stat_info_sub)
        symmetry      = &
             control_get_symmetry(this%ctrl,stat_info_sub)
        dynamic_density_ref = &
             control_get_dynamic_density_ref(this%ctrl,stat_info_sub)
        Newtonian     = &
             control_get_Newtonian(this%ctrl,stat_info_sub)
        Brownian      = &
             control_get_Brownian(this%ctrl,stat_info_sub)    
        stress_tensor = &
             control_get_stress_tensor(this%ctrl,stat_info_sub)
        
#ifdef __PARTICLES_STRESS_SEPARATE
        stress_tensor_p = stress_tensor
        stress_tensor_v = stress_tensor
        stress_tensor_r = stress_tensor .AND. Brownian
#else
        stress_tensor_p = .FALSE.
        stress_tensor_v = .FALSE.
        stress_tensor_r = .FALSE.
#endif
        p_energy      = &
             control_get_p_energy(this%ctrl,stat_info_sub)      
        flow_v_fixed  = &
             control_get_flow_v_fixed(this%ctrl,stat_info_sub)
        adaptive_dt = &
             control_get_adaptive_dt(this%ctrl,stat_info_sub)
         
        
        !----------------------------------------------------
        ! Physics parameters :
        !----------------------------------------------------
        
        num_species    = &
             physics_get_num_species(this%phys,stat_info_sub)
        num_dim        = &
             physics_get_num_dim(this%phys,stat_info_sub)
        dt             = &
             physics_get_dt(this%phys,stat_info_sub)
        step_start     = &
             physics_get_step_start(this%phys,stat_info_sub)
        step_current   = step_start
        step_end       = &
             physics_get_step_end(this%phys,stat_info_sub)
        time_current   = &
             physics_get_time_start(this%phys,stat_info_sub)
        time_end       = &
             physics_get_time_end(this%phys,stat_info_sub)
        eigen_dynamics = &
             physics_get_eigen_dynamics(this%phys,stat_info_sub)
        
        !----------------------------------------------------
        ! Set dt in right hand side obeject.
        ! where dt is used for random force.
        ! Since dt can be a variable during simulation,
        ! (e.g, dt for relax run is different for real run,
        !  or adaptive time step)
        ! it is necessary to update dt in rhs.
        !----------------------------------------------------
        
        CALL particles_get_rhs(this%particles,trhs,stat_info_sub)
        CALL rhs_set_dt(trhs,dt,stat_info_sub)
        
        !----------------------------------------------------
        ! Check if fixed in-let flow velocity is desirded.
        ! If yes, get frequency which is used to 
        ! adjust body force and achieve flow velocity.
        !----------------------------------------------------
        
        IF( flow_v_fixed ) THEN
           
           flow_adjust_freq = &
                physics_get_flow_adjust_freq(this%phys,stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Colloid parameters :
        !
        ! For complex fluids, return a pointer of colloid.
        !
        ! Set colloidal boundary particles velocity according
        ! to colloid translation and rotation speed.
        !
        ! Note that even there is no colloid, set coll_mom=0
        ! for consistency reason.
        !----------------------------------------------------
        
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF ( num_colloid > 0 ) THEN
           
           CALL physics_get_colloid(this%phys,&
                colloids,stat_info_sub)
           
        END IF
        
        ALLOCATE(coll_mom(num_dim))
        coll_mom(:) = 0.0_MK

        !----------------------------------------------------
        ! Boundary parameters :
        !
        ! If there is a shear, update boundary with current
        ! time. i.e., set initial velocity etc.
        !
        ! Walls' velocity in case of oscillating shear;
        ! Sheared length, in case of Lees-Edwards boundary.
        !----------------------------------------------------
        
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        num_sym      = &
             boundary_get_num_sym(tboundary,stat_info_sub)
        num_wall_sym = &
             boundary_get_num_wall_sym(tboundary,stat_info_sub)
        num_wall_solid = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        num_le       = &
             boundary_get_num_le(tboundary,stat_info_sub)
        num_shear    = &
             boundary_get_num_shear(tboundary,stat_info_sub)
        
        IF ( num_shear > 0 ) THEN
           
           CALL boundary_update_boundary(tboundary,&
                time_current,time_current,stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_marching : ", &
                   "boundary updating boundary failed !"
              stat_info = -1           
              GOTO 9999           
           END IF
           
        END IF ! num_shear > 0
        
#ifdef __DEBUG

#ifdef __DEBUG_MARCHING
        
        !----------------------------------------------------
        ! For debugging, output information.
        !----------------------------------------------------
        
        num_part_real  = &
             particles_get_num_part_real(this%particles,stat_info_sub)
        num_part_all  = &
             particles_get_num_part_all(this%particles,stat_info_sub)
        num_part_ghost = &
             particles_get_num_part_ghost(this%particles,stat_info_sub)
        
        CALL particles_get_m(this%particles,t_rho, &
             num_part_all,stat_info_sub)
        CALL debug_write_output(global_debug,rank,"marching", &
             'm_real',0,t_rho,1,num_part_real,stat_info_sub)
        CALL debug_write_output(global_debug,rank,"marching", &
             'm_all',0,t_rho,1,num_part_all,stat_info_sub)
        
        !----------------------------------------------------
        ! Sychronize every process.
        !----------------------------------------------------
        
        CALL MPI_Barrier(comm,stat_info_sub)        
        
#endif
        
#endif 

#ifdef __FLOW_DEVELOPED
        
        !----------------------------------------------------
        ! For reading-external-particles, no need to set
        ! flow to be developed.
        !----------------------------------------------------

        IF ( .NOT. read_external ) THEN
           
           CALL particles_set_flow_developed(this%particles,stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_marching : ", &
                   "Setting particles flow developed failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
           IF ( num_colloid > 0 ) THEN
              
              !CALL colloid_set_flow_developed(colloids,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_marching : ", &
                      "Setting colloid flow developed failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF ! num_colloid
        
        END IF ! read_external
        
#endif        
        !----------------------------------------------------
      	! Create ghost particles on righ-inner(-up)
        ! boundaries (symmetry) or all boundaries 
        ! (non-symmetry) for each process.
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
           PRINT *, "marching_marching : ", &
                "Creating ghosts with x, m, IDs failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! According to different boundary condition,
        ! we have set up the IDs of ghost particles,
        ! which are boundary particles also.
        !----------------------------------------------------
        
        CALL particles_set_boundary_ghost_id(this%particles, &
             stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "marching_marching : ", &
                "Setting boundary ghost ID failed !"
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
      	! to build neighbor list(e.g., cell list).
      	!----------------------------------------------------
        
        CALL particles_get_x(this%particles,t_x, &
             num_part_all,stat_info_sub)        
        CALL technique_build_list(this%tech,t_x, &
             num_part_all,symmetry,stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN
           PRINT *,"marching_marching : ", &
                "Building neighbour list failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! Compute mass density/number density at time_current.
        !----------------------------------------------------
        
        CALL particles_compute_density(this%particles,stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN           
           PRINT *, "marching_marching : ", &
                "Computing density failed !"
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
           
           IF (stat_info_sub /= 0) THEN           
              PRINT *, "marching_marching : ",&
                   "Receiveing contribution from ghosts failed !"
              stat_info = -1
              GOTO 9999           
           END IF
           
        END IF
        
          
        IF ( num_colloid > 0 ) THEN
           
           CALL particles_set_colloid_velocity(this%particles,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'marching_marching : ',&
                   'Setting colloid velocity failed !'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF

        IF ( num_shear > 0 ) THEN
           
           CALL particles_set_boundary_velocity(this%particles,&
                stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_marching : ", &
                   "particles setting boundary failed !"
              stat_info = -1           
              GOTO 9999           
           END IF
           
        END IF ! num_shear > 0
        
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
             l_map_ct =(.NOT. Newtonian), &
             stat_info=stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN
           PRINT *, 'marching_marching : ', &
                'Updating ghosts with x, rho, v failed !'
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! The number of real, ghost and all particles
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
        ! boundary conditions, we have to set up veloicity 
        ! of boundary particles, which ghost particles.
        !----------------------------------------------------
        
        CALL particles_set_boundary_ghost_velocity(this%particles, &
             stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN           
           PRINT *, "marching_marching : ", &
                "Setting boundary ghosts velocity failed !"
           stat_info = -1
           GOTO 9999           
        END IF
        
        
        !----------------------------------------------------
        ! If dynamic rho_ref is required,
        ! we find the minimum density during simulation
        ! set by particles_set_rho_ref into state equation.
        !----------------------------------------------------
        
        IF ( dynamic_density_ref ) THEN
           
           CALL particles_find_density_extreme(this%particles, &
                comm, MPI_PREC, stat_info_sub)
           
           IF(stat_info_sub /=0) THEN
              PRINT *, "marching_marching : ", &
                   "Finding density extrem failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
           CALL particles_set_stateEquation_rho_ref(this%particles, &
                stat_info_sub)
           
           IF(stat_info_sub /=0) THEN
              PRINT *, "marching_marching : ", &
                   "Setting density extreme failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
           rho_min = particles_get_rho_min(this%particles, &
                stat_info_sub)
           rho_max = particles_get_rho_max(this%particles, &
                stat_info_sub)
           
           !-------------------------------------------------
           ! Save rho_min and rho_max as statistic data also.
           !-------------------------------------------------
           
           CALL statistic_set_rho_min(this%statis,rho_min, &
                stat_info_sub)
           
           CALL statistic_set_rho_max(this%statis,rho_max, &
                stat_info_sub)
           
        END IF ! dynamic_density_ref
        
        !----------------------------------------------------
      	! Compute pressure as scalar variable
        ! (including ghost particles) at time_current.
        ! since the pressure of ghost particles 
        ! are used for calculating the interaction(force) also.
        !----------------------------------------------------
        
        CALL particles_compute_pressure(this%particles,&
             num_part_all,stat_info_sub)
        
        IF(stat_info_sub /=0) THEN
           PRINT *, "marching_marching : ", &
                "Computing pressure failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
      	! Compute pressure tensor for non-Newtonian case.
        ! (including ghost particles) at time_current.
        !----------------------------------------------------
        
        IF  ( .NOT. Newtonian ) THEN
           
           CALL particles_compute_pressure_tensor(this%particles,&
                num_part_all,stat_info_sub)
           
           IF(stat_info_sub /=0) THEN
              PRINT *, "marching_marching : ", &
                   "Computing pressure tensor failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
      	! Compute any interaction between particles 
        ! at time_current.
        ! e.g., force, velocity gradient tensor, 
        ! stress tensor etc.
      	!----------------------------------------------------
        
        CALL particles_compute_interaction(this%particles, &
             stat_info_sub)
        
        IF(stat_info_sub /=0 ) THEN
           PRINT *, "marching_marching : ", &
                "Computing interaction failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! For symmetry inter-process communication:
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
           
#ifdef __PARTICLES_FORCE_SEPARATE
           CALL particles_map_ghost_put(this%particles, &
                l_map_x   = .TRUE., &
                l_map_f   = .TRUE.,  &
                l_map_fp  = .TRUE., &
                l_map_fv  = .TRUE.,  &
                l_map_fr  = Brownian,  &
                l_map_s   = stress_tensor,   &
                l_map_sp  = stress_tensor_p, &
                l_map_sv  = stress_tensor_v, &
                l_map_sr  = stress_tensor_r, &
                l_map_vgt = (.NOT. Newtonian),&
                l_map_au  = p_energy, &
                stat_info = stat_info_sub)
           
#else
           CALL particles_map_ghost_put(this%particles, &
                l_map_x   = .TRUE., l_map_f  = .TRUE., &
                l_map_s   = stress_tensor,    &
                l_map_vgt = (.NOT. Newtonian),&
                l_map_au  = p_energy, &
                stat_info = stat_info_sub)
#endif

           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, 'marching_marching : ',&
                   'Receiving force(stress, vgt,au) from ghosts failed !'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! symmetry
        
        
        !----------------------------------------------------
        ! If the in-let flow velocity is a input parameter,
        ! we adjust body force of fluid(and colloid)
        ! to achieve desired flow velocity.
        !----------------------------------------------------
        
        CALL physics_get_body_force(this%phys, &
             body_force,stat_info_sub)
        
        IF( flow_v_fixed ) THEN
           
           CALL physics_get_body_force_d(this%phys, &
                body_force_d,stat_info_sub)
           
           body_force(1:num_dim) = body_force(1:num_dim) + &
                body_force_d(1:num_dim)
           
           CALL physics_set_body_force(this%phys, &
                body_force,stat_info_sub)
           
           IF ( num_colloid > 0 ) THEN
              
              CALL colloid_set_body_force(colloids, &
                   body_force,stat_info_sub)
              
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Apply body forces to real fluid particles.
        !----------------------------------------------------
        
        CALL particles_apply_body_force(this%particles,&
             num_part_real,stat_info_sub)
        
        IF(stat_info_sub /=0 ) THEN
           PRINT *, "marching_marching : ",&
                "Applying body force failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! In case of Non-Newtonian fluid,
        ! we use Oldroyd-B model.
        !----------------------------------------------------
        
        IF( .NOT. Newtonian ) THEN
           
           !-------------------------------------------------
           ! In case of egenvector dynamics,
           ! compute accleration of eigenvalues and
           ! eigenvectors, then integrate them.
           ! Finally compute conformation tensor.
           !-------------------------------------------------
           
           IF ( eigen_dynamics ) THEN
              
              !----------------------------------------------
              ! Compute matrix element of eigen-dynamics 
              ! from velocity gradient tensor.
              !----------------------------------------------
              
              CALL particles_compute_evgt(this%particles,&
                   num_part_real,stat_info_sub)
              
              !----------------------------------------------
              ! Compute accelerations of eigenvalues.
              !----------------------------------------------
              
              CALL particles_compute_aeval(this%particles,&
                   num_part_real,stat_info_sub)
              
              IF (stat_info_sub /= 0 ) THEN
                 PRINT *, &
                      "marching_marching : ",&
                      "Computing aeval failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Compute accelerations of eigenvectors.
              !----------------------------------------------
              
              CALL particles_compute_aevec(this%particles,&
                   num_part_real,stat_info_sub)
              
              IF (stat_info_sub /= 0 ) THEN
                 PRINT *,&
                      "marching_marching : ",&
                      "Computing aevec failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           ELSE ! evolution of conformation tensor.

              !----------------------------------------------
              ! Compute accleration of conformation
              ! tensor from velocity gradient tensor.
              !----------------------------------------------
              
              CALL particles_compute_act(this%particles,&
                   num_part_real,stat_info_sub)
              
              IF (stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_marching : ",&
                      "Computing act failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF ! eigen-dynamics
           
        END IF  ! non-Newtonian
        
        
        !----------------------------------------------------
        ! Since paire-wise force have been calculated,
        ! sum up the force on colloids,
        ! if colloids are present.
        ! Compute colloid-colloid and colloid-wall 
        ! interactions.
        !----------------------------------------------------
        
        IF ( num_colloid > 0 ) THEN
           
           !-------------------------------------------------
           ! Sum up force/torque exerted on boundary particles
           ! of colloids on this process.
           !-------------------------------------------------
           
           CALL particles_collect_colloid_interaction(this%particles,&
                coll_drag,coll_torque,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_marching : ",&
                   "Summing up interaction on colloid locally has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Collecting force/torque on colloids by 
           ! boundary particles from all processes.
           !-------------------------------------------------
           
           CALL colloid_collect_particles_interaction(colloids,&
                comm,MPI_PREC,coll_drag,coll_torque,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_marching : ",&
                   "Summing up interaction on colloid globally has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Compute colloid-colloid and colloid-wall
           ! interactions, return drag on the wall.
           !-------------------------------------------------
           
           CALL colloid_compute_interaction(colloids,comm,&
                MPI_PREC,wall_drag_c(1:num_dim,1:num_dim*2),&
                stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_marching : ",&
                   "c-c or c-w interaction has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Apply body force on colloids.
           !-------------------------------------------------
           
           CALL colloid_apply_body_force(colloids,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_marching : ",&
                   "Applying body force on colloids has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Compute colloids accelerations, i.e.,
           ! translation and rotation.
           !-------------------------------------------------
           
           CALL colloid_compute_acceleration(colloids,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_marching : ",&
                   "Computing colloids accelerations has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Compute colloids' statistics first, if present.
           ! Compute initial statistics of all real particles.
           !-------------------------------------------------
           
           CALL colloid_compute_statistic(colloids,stat_info_sub)
           coll_k = colloid_get_k_energy_tot(colloids,stat_info_sub)
           CALL colloid_get_mom_tot(colloids,coll_mom,stat_info_sub)
           
        END If ! num_colloid > 0
        
        
        CALL statistic_compute_statistic(this%statis, &
             this%particles,coll_k,coll_mom,stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN
           PRINT *, "marching_marching : ",&
                "Computing statistics failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check if we need to adapt body force to get
        ! desired flow velocity.
        !----------------------------------------------------
        
        IF( flow_v_fixed ) THEN
           
           CALL marching_adjust_flow_v(this,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN              
              PRINT *, "marching_marching : ",&
                   "Adjusting flow velocity failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! If there is wall using symmetry or solid boundary
        ! particles, or Lees-Edwards boundary particles.
        ! we sum up the interaction on boundaries.
        !----------------------------------------------------
        
        IF ( num_shear > 0 ) THEN
           
#ifdef __WALL_FORCE_SEPARATE
           
           CALL particles_collect_boundary_interaction(this%particles,&
                wall_drag_p(1:num_dim,1:num_dim*2),& 
                wall_drag_pp(1:num_dim,1:num_dim*2),& 
                wall_drag_pv(1:num_dim,1:num_dim*2),& 
                wall_drag_pr(1:num_dim,1:num_dim*2),& 
                stat_info_sub)

#else
           
           CALL particles_collect_boundary_interaction(this%particles,&
                wall_drag_p(1:num_dim,1:num_dim*2), stat_info=stat_info_sub)

#endif
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_marching : ",&
                   "Summing up interaction on boundary locally has problem !"
              stat_info = -1
              GOTO 9999
           END IF
          
#ifdef __WALL_FORCE_SEPARATE

           CAll boundary_collect_particles_interaction(tboundary,comm,&
                MPI_PREC, wall_drag_p(1:num_dim,1:num_dim*2),&
                wall_drag_pp(1:num_dim,1:num_dim*2),& 
                wall_drag_pv(1:num_dim,1:num_dim*2),& 
                wall_drag_pr(1:num_dim,1:num_dim*2),&             
                stat_info_sub)
#else 
           CAll boundary_collect_particles_interaction(tboundary,comm,&
                MPI_PREC, wall_drag_p(1:num_dim,1:num_dim*2),&
                stat_info=stat_info_sub)
#endif
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_marching : ",&
                   "Summing up particles contribution on boundary has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           CAll boundary_collect_colloid_interaction(tboundary,comm,&
                MPI_PREC, wall_drag_c(1:num_dim,1:num_dim*2),stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_marching : ",&
                   "Summing up colloids contribution on boundary has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Reset boundary particles interaction,
           ! both real and ghost particles.
           !-------------------------------------------------
           
           CALL particles_reset_boundary_interaction(this%particles, &
                stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_marching : ",&
                   "Resetting boundary particles interaction failed !"
              stat_info = -1
              GOTO 9999
           END IF

           CALL particles_reset_boundary_ghost_interaction(this%particles, &
                stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_marching : ",&
                   "Resetting boundary ghost particles interaction failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! num_shear > 0
        
        
        IF ( rank == 0 ) THEN
           
           !-------------------------------------------------
           ! Open files, e.g., statistics, boundary files.
           !-------------------------------------------------
           
           CALL io_open(this%io,rank,&
                1,num_shear,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_marching : ",&
                   "Opening (s,b) files failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
	! Write initial conditions into output files
      	!----------------------------------------------------
        
        CALL io_write_condition_check(this%io,&
             step_current,time_current,&
             write_particles=l_write_particles,stat_info=stat_info_sub)
        
        num_part_real = &
             particles_get_num_part_real(this%particles,stat_info_sub)
        
        CALL io_write(this%io,rank,step_current,time_current,&
             this%particles,num_part_real,&
             this%statis,stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN              
           PRINT *, "marching_marching : ",&
                "Initial writing failed !"
           stat_info = -1
           GOTO 9999
        END IF
        

#ifdef __DEBUG

#ifdef __DEBUG_INTEGRATE
        
        !----------------------------------------------------
        ! For debug purpose, open file to record time cost.
        !----------------------------------------------------
        
        IF ( rank == 0 .AND. debug_flag == 3) THEN
           
           CALL debug_open_time(global_debug,rank,stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_marching : ", &
                   "debug_open_time has problem ! "
              stat_info_sub = -1
              GOTO 9999
           END IF
           
        END IF
        
#endif

#endif        
        !----------------------------------------------------
        !----------------------------------------------------
        ! ##################################################
        ! ######      Time loop for integration.      ######
        ! ##################################################
        !----------------------------------------------------
      	!----------------------------------------------------

        DO WHILE ( (time_end >=0.0_MK .AND. &
             time_current <  time_end ) .OR. &
             ( step_end >=0 .AND. &
             step_current <  step_end ) )
           
        
           !-------------------------------------------------
           ! For adaptive dt, calculate dt now and update it
           !-------------------------------------------------
           
           IF ( adaptive_dt > 0 ) THEN
              
              !----------------------------------------------
              ! Find the maximum acceleration of 
              ! fluid particle
              !----------------------------------------------
              
              CALL particles_find_force_extreme(this%particles,&
                   comm, MPI_PREC,stat_info_sub)
              CALL particles_compute_dt_f(this%particles,stat_info_sub)
              dt_f = &
                   particles_get_dt_f(this%particles,stat_info_sub)
              
              IF ( num_colloid > 0 ) THEN
                 
                 !-------------------------------------------
                 ! Find the maximum acceleration of 
                 ! colloidal particle.
                 !-------------------------------------------
            
                 CALL colloid_find_force_extreme(colloids,stat_info_sub)
                 CALL colloid_compute_dt_f(colloids,stat_info_sub)
                 
                 dt_f_c = &
                      colloid_get_dt_f(colloids,stat_info_sub)
                 
                 IF ( dt_f_c > 0.0_MK ) THEN
                    
                    IF ( dt_f <= 0.0_MK ) THEN
                       
                       dt_f = dt_f_c
                       
                    ELSE IF ( dt_f_c < dt_f ) THEN
                       
                       dt_f = dt_f_c
                       
                    END IF
                    
                 END IF
                 
              END IF
              
              CALL physics_adapt_dt(this%phys,dt_f,stat_info_sub)
              
              dt  = &
                   physics_get_dt(this%phys,stat_info_sub)
              CALL rhs_set_dt(trhs,dt,stat_info_sub)
              
           END IF
           
           !-------------------------------------------------
           ! Integrate with time
           !-------------------------------------------------
           
           CALL marching_integrate(this,time_current,dt,&
                stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN              
              PRINT *, "marching_marching : ",&
                   "Integrating failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Compute new images(position and velocity)
           ! of colloids.
           !-------------------------------------------------
           
           IF ( num_colloid > 0 ) THEN
              
              CALL colloid_compute_image(colloids,stat_info_sub)
              
              IF( stat_info_sub /=0 ) THEN              
                 PRINT *, "marching_marching : ",&
                      "colloid computing image failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF
           
           !-------------------------------------------------
           ! Update current step and time.
           !-------------------------------------------------
           
           step_current = step_current + 1
           time_current = time_current + dt
           
           CALL physics_set_step_current(this%phys,step_current,&
                stat_info_sub)
           CALL physics_set_time_current(this%phys,time_current,&
                stat_info_sub)
          
           !-------------------------------------------------
           ! Check what we need to write at this time step.
           !-------------------------------------------------
           
           CALL io_write_condition_check(this%io,&
                step_current,time_current,&
                write_statistic=l_write_statistic,&
                stat_info=stat_info_sub)
           
           IF(stat_info_sub /=0 ) THEN              
              PRINT *,"marching_marching : ",&
                   "Checking writing conditions failed!"              
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! If we need to write statistic
           ! information, we have to 
           ! calculate first.
           !-------------------------------------------------
           
           IF( l_write_statistic ) THEN
              
              !----------------------------------------------
              ! Compute colloids' statistics first,
              ! if present.           
              ! Compute statistics of real particles.
              !----------------------------------------------
              
              IF( num_colloid > 0) THEN
                 
                 CALL colloid_compute_statistic(colloids,stat_info_sub)
                 coll_k = colloid_get_k_energy_tot(colloids,stat_info_sub)
                 CALL colloid_get_mom_tot(colloids,coll_mom,stat_info_sub)
                 
              END IF
              
              CALL statistic_compute_statistic(this%statis, &
                   this%particles,coll_k,coll_mom,stat_info_sub)
              
              IF(stat_info_sub /=0 ) THEN
                 PRINT *, "marching_marching : ", &
                      "Computing statistics failed !"
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
              
           END IF
           
           !-------------------------------------------------
           ! Write informations into output files
           !-------------------------------------------------
           
           num_part_real = &
                particles_get_num_part_real(this%particles,stat_info_sub)
           
           CALL io_write(this%io,rank,step_current,time_current,&
                this%particles,num_part_real,&
                this%statis,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN              
              PRINT *, "marching_marching : ",&
                   "Step writing failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Check if we need to adapt body force to get
           ! desired flow velocity.
           !-------------------------------------------------
           
           IF( flow_v_fixed .AND. &
                MOD(step_current,flow_adjust_freq)==0 ) THEN
              
              CALL marching_adjust_flow_v(this,stat_info_sub)
              
              IF( stat_info_sub /=0 ) THEN              
                 PRINT *, "marching_marching : ",&
                      "Adjusting flow velocity failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF
           
        END DO ! time_current < time_end
        
        !----------------------------------------------------
        !----------------------------------------------------
        !####################################################
        !#######   End of time loop for integration.  #######
        !####################################################
        !----------------------------------------------------
        !----------------------------------------------------
        

#ifdef __DEBUG
        
        IF ( rank == 0 .AND. debug_flag == 3) THEN
           
           CALL debug_close_time(global_debug,rank,stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_marching : ", &
                   "debug_close_time has problem ! "
              stat_info_sub = -1
              GOTO 9999
           END IF

        END IF
        
#endif
        
        !----------------------------------------------------
        ! Writing everything in the end is always useful.
        !
        ! First check if the final step has already 
        ! been written during last step of the loop.
        !
        ! If not, output the final step with 
        ! all information, such as particles configuration,
        ! conformation tensor, statistic file, colloid file, 
        ! restart physics config, 
        ! restart particles configuration,
        ! restart conformation tensor.
        ! We do so by flipping over the write_condition flags.
        !----------------------------------------------------
        
        CALL io_write_condition_check(this%io,&
             step_current,time_current,&
             l_write_particles,l_write_conformation,&
             l_write_colloid,&
             l_write_statistic,l_write_boundary,&
             l_write_restart_physics, &
             l_write_restart_particles, &
             l_write_restart_conformation, &
             stat_info_sub)
        
        IF( .NOT. l_write_statistic ) THEN
           
           !-------------------------------------------------
           ! Compute colloids' statistics first, if present.
           !-------------------------------------------------
           
           IF( num_colloid > 0) THEN
              
              CALL colloid_compute_statistic(colloids,stat_info_sub)
              coll_k = colloid_get_k_energy_tot(colloids,stat_info_sub)
              CALL colloid_get_mom_tot(colloids,coll_mom,stat_info_sub)
              
           END IF
           
           CALL statistic_compute_statistic(this%statis, &
                this%particles,coll_k,coll_mom,stat_info_sub)
           
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
              
        END IF
        
        CALL io_write_condition_set(this%io,&
             .NOT. l_write_particles, &
             .NOT. l_write_conformation,&
             .NOT. l_write_colloid,&
             .NOT. l_write_statistic,&
             .NOT. l_write_boundary,&
             .NOT. l_write_restart_physics  ,&
             .NOT. l_write_restart_particles ,&
             .NOT. l_write_restart_conformation,&
             stat_info_sub)
        
        !----------------------------------------------------
        ! Write final step output.
        !----------------------------------------------------
        
        num_part_real = &
             particles_get_num_part_real(this%particles,stat_info_sub)
        
        CALL io_write(this%io,rank,step_current,time_current,&
             this%particles,num_part_real,&
             this%statis,stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN              
           PRINT *, "marching_marching : ",&
                "Final writing failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Close statistic,boundary and colloid files.
        !----------------------------------------------------
        
        CALL io_close(this%io,1,num_shear,stat_info_sub)

        IF( stat_info_sub /=0 ) THEN           
           PRINT *,"marching_marching : ",&
                "Closing (s,b,c) files failed !"
           stat_info = -1           
        END IF
        
        !----------------------------------------------------
        ! Release all dynamics memories.
        !----------------------------------------------------
        
        IF(ASSOCIATED(body_force)) THEN
           DEALLOCATE(body_force)
        END IF

        IF(ASSOCIATED(body_force_d)) THEN
           DEALLOCATE(body_force_d)
        END IF
        
        IF(ASSOCIATED(coll_mom)) THEN
           DEALLOCATE(coll_mom)
        END IF
        
        IF(ASSOCIATED(coll_drag)) THEN
           DEALLOCATE(coll_drag)
        END IF

        IF(ASSOCIATED(coll_torque)) THEN
           DEALLOCATE(coll_torque)
        END IF
        
        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF
      
        IF(ASSOCIATED(t_x)) THEN
           DEALLOCATE(t_x)
        END IF
        
        
#ifdef __DEBUG
        
        IF( debug_flag == 2) THEN
           CALL debug_substop(global_debug,rank, &
                'marching_marching :', &
                time_routine_start, stat_info_sub)
        END IF 

        IF(ASSOCIATED(t_rho)) THEN
           DEALLOCATE(t_rho)
        END IF

        IF(ASSOCIATED(t_f)) THEN
           DEALLOCATE(t_f)
        END IF
        
        IF(ASSOCIATED(t_id)) THEN
           DEALLOCATE(t_id)
        END IF
        
#endif
        
        RETURN
        
        
      END SUBROUTINE marching_marching
      
      
      
