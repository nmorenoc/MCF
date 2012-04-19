      SUBROUTINE marching_integrate_VV(this,step,time,dt,stat_info)
        !----------------------------------------------------
        ! Subroutine  : marching_integrate_VV
        !----------------------------------------------------
        !
        ! Purpose     : For time integration using Velocity
        !               Verlet with symmetric/non-symmetry 
        !               inter-processes communication.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.4 30.07 2009,
        !               include non-Newtonian viscoelastic
        !               Oldroyd-B mode.
        !
        !               V0.3 23.07 2009, merge
        !               marching_integer_vv_nonsym() and
        !               marching_integer_vv_sym() together
        !               
        !               V0.2 09.07 2009, 
        !               check again the work flow and
        !               supply with more comments.
        !
        !               V0.1 30.06 2009, original version.
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
        ! time       : current time.
        ! dt         : time step.
        ! stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Marching), INTENT(INOUT)   :: this
        INTEGER, INTENT(IN)             :: step
        REAL(MK), INTENT(IN)            :: time
        REAL(MK), INTENT(IN)            :: dt
        INTEGER, INTENT(OUT)	        :: stat_info

	!----------------------------------------------------
	! Local variables start here :
	!----------------------------------------------------
	
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Control parameters :
     	!----------------------------------------------------
        
        LOGICAL                         :: dynamic_density_ref
        LOGICAL                         :: symmetry
        LOGICAL                         :: Newtonian
        LOGICAL                         :: Brownian
        LOGICAL                         :: stress_tensor
        LOGICAL                         :: stress_tensor_p
        LOGICAL                         :: stress_tensor_v
        LOGICAL                         :: stress_tensor_r     
        LOGICAL                         :: p_energy
        
        !----------------------------------------------------
        ! Physics parameters.(colloids)
     	!----------------------------------------------------
        
        INTEGER                         :: num_species
        INTEGER                         :: num_dim, i
        INTEGER                         :: step_start

        LOGICAL                         :: eigen_dynamics
        
        INTEGER                         :: num_colloid
        INTEGER                         :: coll_sub_time_step
        REAL(MK)                        :: dt_sub_time_step        
        TYPE(Colloid), POINTER          :: colloids
        REAL(MK), ALLOCATABLE, DIMENSION(:,:)   :: coll_drag
        REAL(MK), ALLOCATABLE, DIMENSION(:,:)   :: coll_torque
        
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: num_wall_solid
        INTEGER                         :: num_shear
        REAL(MK),DIMENSION(3,6)         :: wall_drag_p
        REAL(MK),DIMENSION(3,6)         :: wall_drag_c
#ifdef __WALL_FORCE_SEPARATE
        REAL(MK),DIMENSION(3,6)         :: wall_drag_pp
        REAL(MK),DIMENSION(3,6)         :: wall_drag_pv
        REAL(MK),DIMENSION(3,6)         :: wall_drag_pr        
#endif
        !----------------------------------------------------
        ! t_x  : position of particles
        !----------------------------------------------------
        
        REAL(MK), DIMENSION(:,:), POINTER       :: t_x 
        
        !----------------------------------------------------
        ! Number of real, all and ghost particles.
     	!----------------------------------------------------
        
        INTEGER                         :: num_part_real
        INTEGER                         :: num_part_all
        INTEGER                         :: num_part_ghost
        
        !----------------------------------------------------
        !  MPI parameters.
        !----------------------------------------------------
        
        INTEGER                         :: rank
        INTEGER                         :: comm
        INTEGER                         :: MPI_PREC
        
#ifdef __DEBUG_INTEGRATE_VV
        
        !----------------------------------------------------
        !  Debug variables
        !----------------------------------------------------
        INTEGER                         :: debug_flag
        REAL(MK)			:: debug_time_start
        REAL(MK)                        :: debug_time0
        REAL(MK)                        :: debug_time1
        REAL(MK), DIMENSION(15)         :: debug_time_record
        INTEGER                         :: debug_index
        
#endif 
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        
        num_colloid = 0
        NULLIFY(colloids)
         
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
        ! Control parameters.
        !----------------------------------------------------
        
        dynamic_density_ref = &
             control_get_dynamic_density_ref(this%ctrl,stat_info_sub)
        symmetry  = &
             control_get_symmetry(this%ctrl,stat_info_sub)
        Newtonian = &
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
        p_energy  = &
             control_get_p_energy(this%ctrl,stat_info_sub) 
        
        !----------------------------------------------------
        ! Physics parameters.
        !----------------------------------------------------
        
        num_species    = &
             physics_get_num_species(this%phys,stat_info_sub)
        num_dim        = &
             physics_get_num_dim(this%phys,stat_info_sub)
        step_start     = &
             physics_get_step_start(this%phys,stat_info_sub)
        eigen_dynamics = &
             physics_get_eigen_dynamics(this%phys,stat_info_sub)
        
        !----------------------------------------------------
        ! Return the object pointer of Class Colloid.
        !----------------------------------------------------
        
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF ( num_colloid > 0 ) THEN           
           
           CALL physics_get_colloid(this%phys,colloids,stat_info_sub)
        
           coll_sub_time_step = &
                colloid_get_sub_time_step(colloids,stat_info_sub)
           dt_sub_time_step   = dt / coll_sub_time_step
           
            ALLOCATE(coll_drag(num_dim,num_colloid))
            ALLOCATE(coll_torque(3,num_colloid))
           
        END IF
        
        !----------------------------------------------------
        ! Get boundary conditions.
        !----------------------------------------------------
        
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        num_wall_solid = & 
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        num_shear      = &
             boundary_get_num_shear(tboundary,stat_info_sub)
    
        !----------------------------------------------------
        ! Number of real, all and ghost particles.
        !----------------------------------------------------
     
        num_part_real  = &
             particles_get_num_part_real(this%particles,stat_info_sub)
        num_part_ghost = &
             particles_get_num_part_ghost(this%particles,stat_info_sub)
        num_part_all   = &
             particles_get_num_part_all(this%particles,stat_info_sub)
        
        !----------------------------------------------------
        ! MPI parameters.
        !----------------------------------------------------
        
        rank     = technique_get_rank(this%tech,stat_info_sub)
        MPI_PREC = technique_get_MPI_PREC(this%tech,stat_info_sub)
        comm     = technique_get_comm(this%tech,stat_info_sub)
        
        
#ifdef __DEBUG_INTEGRATE_VV
        debug_flag = &
             debug_get_flag(global_debug,stat_info_sub)
        SELECT CASE (debug_flag)
        CASE (2)
           CALL debug_substart(global_debug,rank, &
                "marching_integrate_VV", &
                debug_time_start,stat_info_sub)
        CASE (3)
           debug_time_start = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = 1
        END SELECT
#endif
        
        !----------------------------------------------------
	! Calculate new velocity 
        ! v'(t+dt) = v(t) + lamda * f(t) * dt,
        ! where lambda=0.5 here.
        !
        ! Remark : In principle, after updating positions,
        !          mapping ghosts should be done immediately,
        !          since ghosts needed to be refreshed.
        !          However, if updating velocity first, 
        !          before mapping ghost, then we don't need
        !          to map force for ghosts, it saves
        !          expense of communication.
        !          The velocities have anyway to be mapped
        !          for ghosts.
	!----------------------------------------------------
        
        
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif
        
        CALL particles_integrate_velocity(this%particles, &
             num_part_real,dt,0.5_MK,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "marching_integrate_VV : ", &
                "Integrating velocity failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
	!----------------------------------------------------
	! Position integration r(t+dt) = 
        ! r(t) + v(t) * dt  + 0.5*f(t)*dt**2,
        ! howver, since we integrated v first already above
        ! as v'(t+dt) = v(t) + 0.5 * f(t) * dt,
        ! here we should have r(t+dt) = r(t) + v'(t)*dt
	!----------------------------------------------------
        
        CALL particles_integrate_position(this%particles, &
             num_part_real,dt,1,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "marching_integrate_VV : ", &
                "Integrating position failed ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif        
        
        !----------------------------------------------------
        ! For non-Newtonian viscoelastic Oldroyd-B model.
        !----------------------------------------------------     
        
        IF( .NOT. Newtonian ) THEN
           
           !-------------------------------------------------
           ! In case of eigen-dynamics, use accelerations
           ! of eigenvalues and eigenvectors to integrate
           ! eigenvalues and eigenvectors.
           ! Finally compute conformation tensor using
           ! evals and evecs.
           !-------------------------------------------------
           
           IF ( eigen_dynamics ) THEN
              
              !----------------------------------------------
              ! Integrate eigenvalues with same order
              ! as velocity, i.e., lamda = 0.5.
              !----------------------------------------------
              
              CALL particles_integrate_eval(this%particles, &
                   num_part_real,dt,0.5_MK,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *,&
                      "marching_integrate_VV : ", &
                      "Integrating eval failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Integrate eigenvectors with same order
              ! as velocity, i.e., lamda = 0.5.
              !----------------------------------------------
              
              CALL particles_integrate_evec(this%particles, &
                   num_part_real,dt,0.5_MK,stat_info_sub)
              
              IF (stat_info_sub /= 0 ) THEN
                 PRINT *,&
                      "marching_integrate_VV : ", &
                      "Integrating evec failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Compute conformation tensor out of
              ! eigenvalues and eigenvectors.
              !----------------------------------------------
              
              CALL particles_compute_ct(this%particles, &
                   num_part_real,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV : ", &
                      "Computing ct failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           ELSE ! evolution of conformation tensor.
              
              !----------------------------------------------
              ! Use acceleration of conformation tensor to
              ! integrate conformation tensor with same
              ! order as velocity, i.e., lamda=0.5.
              !----------------------------------------------
              
              CALL particles_integrate_ct(this%particles, &
                   num_part_real,dt,0.5_MK,stat_info_sub)
              
              IF (stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV: ", &
                      "Integrating ct failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF ! eigen_dynamics
           
        END IF ! non-Newtonian
        
        
        !----------------------------------------------------
        ! Check if potential energy is needed.
        !----------------------------------------------------
        
        IF( p_energy ) THEN
           
           CALL particles_integrate_potential_energy(this%particles, &
                num_part_real,dt,0.5_MK,stat_info_sub)
           
           IF (stat_info_sub /= 0 ) THEN
              PRINT *, "marching_integrate_VV: ", &
                   "Integrating potential energy failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! p_energy
        
        
       
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif

        IF ( num_colloid > 0 ) THEN           
           

           !-------------------------------------------------
           ! Sum up the force of solvent particles on 
           ! parts of colloids from local processor.
           !-------------------------------------------------
           
           CALL particles_collect_colloid_interaction(&
                this%particles,coll_drag,coll_torque,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_integrate_VV: ",&
                   "Summing up interaction on colloid locally has problem!"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Sum up force/torque of solvent partilces
           ! exerted on colloids from all processes.
           !-------------------------------------------------
           
           CALL colloid_collect_particles_interaction(colloids,&
                comm,MPI_PREC,coll_drag,coll_torque,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_integrate_VV: ",&
                   "Summing up interaction on colloid globally has problem!"
              stat_info = -1
              GOTO 9999
           END IF
           
           
           DO i = 1, coll_sub_time_step
              
              !----------------------------------------------
              ! Compute the rotation vector using rotating
              ! velocity, using desired order.
              !----------------------------------------------
        
              CALL colloid_compute_rotation_vector(colloids,&
                   step-1+i-step_start,dt_sub_time_step,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV: ", &
                      "Computing rotation vector of colloids failed!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Compute rotation matrix from rotation vector.
              !----------------------------------------------
           
              CALL colloid_compute_rotation_matrix(colloids,stat_info_sub)
              CALL colloid_compute_accumulation_matrix(colloids,stat_info_sub)
              CALL colloid_compute_accumulation_vector(colloids,stat_info_sub)
              
              IF ( stat_info_sub /=0 ) THEN
                 
                 PRINT *, "marching_integrate_VV: ", &
                      "Computing rotaiton matrix failed! "
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Compute colloid boundary particle's new relative
              ! position to the colloid center after rotation.
              !----------------------------------------------
              
              CALL particles_compute_colloid_relative_position(&
                   this%particles,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV: ", &
                      "Computing boundary particles relative position failed!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Integrate the positions of all colloids'
              ! centers using desired order.
              !----------------------------------------------
              
              CALL colloid_integrate_position(colloids,&
                   step-1+i-step_start,dt_sub_time_step,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV: ", &
                      "Integrating colloids position failed!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              
              !----------------------------------------------
              ! Compute colloid boundary particle's new
              ! absolute position after the colloid center
              ! is updated.
              !----------------------------------------------
              
              CALL particles_compute_colloid_absolute_position(&
                   this%particles,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV: ", &
                      "Computing boundary particles absolute position failed!"
                 stat_info = -1
                 GOTO 9999
              END IF
           
              !-------------------------------------------------
              ! Integrate velocity using desired accuracy order.
              !-------------------------------------------------
              
              CALL colloid_integrate_velocity(colloids,&
                   step-1+i-step_start,dt_sub_time_step,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV: ", &
                      "Integrating colloids velocity failed!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Add up force/torque from colloid-colloid and
              ! colloid-wall interactions.
              !----------------------------------------------
              
              CALL colloid_compute_interaction(colloids,comm, &
                   MPI_PREC,coll_drag,coll_torque, &
                   wall_drag_c(1:num_dim,1:num_dim*2),stat_info_sub)
              
              IF( stat_info_sub /=0 ) THEN
                 PRINT *, "marching_integrate_VV: ",&
                      "c-c or c-w interaction has problem!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Apply body force on colloids.
              !----------------------------------------------
              
              CALL colloid_apply_body_force(colloids,stat_info_sub)
              
              IF( stat_info_sub /=0 ) THEN
                 PRINT *, "marching_integrate_VV: ", &
                      "Applying body force on colloids has problem!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Compute colloids accelerations, i.e.,
              ! translation and rotation.
              !----------------------------------------------
           
              CALL colloid_compute_acceleration(colloids,stat_info_sub)
              
              IF( stat_info_sub /=0 ) THEN
                 PRINT *, "marching_integrate_VV: ",&
                      "Computing colloids accelerations has problem!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END DO ! i = 1, coll_sub_time_step
           
           !-------------------------------------------------
           ! In case colloids centers go out of physical
           ! boundary, adjust them according to boundary
           ! condition.
           !-------------------------------------------------

           CALL colloid_adjust_colloid(colloids,stat_info_sub)

           IF (stat_info_sub /= 0 ) THEN
              PRINT *, "marching_integrate_VV: ", &
                   "Adjusting colloids r or v failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! num_colloid > 0

#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif   
        
        !----------------------------------------------------
        ! Update boundary :
        !
        ! Walls' velocity in case of oscillating shear;
        ! Sheared length, in case of Lees-Edwards boundary.
        !----------------------------------------------------

        IF ( num_shear > 0 ) THEN
           
           CALL boundary_update_boundary(tboundary, &
                time+dt,time+0.5_MK*dt,stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_integrate_VV : ", &
                   "Updating boundary failed  !"
              stat_info = -1
              GOTO 9999
           END IF
           
           CALL particles_integrate_boundary_position(this%particles, &
                dt,1,stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_integrate_VV : ", &
                   "Integrating boundary position failed ! "
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! num_shear > 0
        
        !----------------------------------------------------
        ! Adjust real particles' r/v after motion,
        ! according to boundary conditions,
        ! in case they go out of the physical domain.
        !----------------------------------------------------

        CALL particles_adjust_particles(this%particles,&
             num_part_real,stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "marching_integrate_VV : ", &
                "Adjusting r or v failed ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! Decompose partially, since positions have changed.
        ! Aussuming the particles can at furthest move to 
        ! neigboring processes at the end of each time step,
        ! we call particles_decompose_partial(), which will
        ! only communicate with neigboring processes.
        !
        ! Density doesn't need to commmunicate, since
        ! it is anyway calculated again.
        ! The same for force ,vgt (velocity gradient tensor)
        ! and potential energy accelearation which will
        ! be allocated and calculated again every step.
        !
        ! Potential energy whether comumincate denpends
        ! on control parameter.
        !----------------------------------------------------

#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif

        CALL particles_decompose_partial(this%particles,&
             l_map_x    = .TRUE., l_map_v  = .TRUE., &
             l_map_m    = .TRUE., l_map_id = .TRUE., &
             l_map_eval = ((.NOT. Newtonian) .AND. eigen_dynamics), &
             l_map_evec = ((.NOT. Newtonian) .AND. eigen_dynamics), &
             l_map_ct   = (.NOT. Newtonian), &
             l_map_u    = p_energy, &
             stat_info  = stat_info_sub)

#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif   

        IF ( stat_info_sub /= 0 ) THEN
           PRINT *,"marching_integrate_VV : ", &
                "Decomposing partially failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! After decompostion, number of real particles
        ! on each process may have changed.
        !----------------------------------------------------
        
        num_part_real = &
             particles_get_num_part_real(this%particles,stat_info_sub)
        num_part_ghost  = &
             particles_get_num_part_ghost(this%particles,stat_info_sub)
        num_part_all  = &
             particles_get_num_part_all(this%particles,stat_info_sub)
      
        
        !----------------------------------------------------
        ! After decomposition of particles, we must create 
        ! new ghosts layers for neighboring processes.
        !
        ! Even using single process, this has to be done,
        ! since it guarantees the boundary condition.
        !  
        ! Has to be done BEFORE building neigbor list.
	!----------------------------------------------------

#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif   

        CALL particles_map_ghost_get(this%particles, &
             l_map_x  = .TRUE., l_map_m =.TRUE., &
             l_map_id = .TRUE., stat_info=stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *,"marching_integrate_VV : ", &
                "Creating ghosts failed!"
           stat_info = -1     
           GOTO 9999
        END IF
        
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif   

        !----------------------------------------------------
        ! According to different boundary condition,
        ! we have set up the IDs of boundary particles
        ! which are ghosts also.
        !----------------------------------------------------
        
        CALL particles_set_boundary_ghost_id(this%particles, &
             stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "marching_integrate_VV: ", &
                "Setting boundary ghosts ID failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! After mapping ghosts number of all(ghosts)
        ! particles on each process might have changed.
        !----------------------------------------------------
        
        num_part_real  = &
             particles_get_num_part_real(this%particles,stat_info_sub)
        num_part_ghost = &
             particles_get_num_part_ghost(this%particles,stat_info_sub)
        num_part_all   = &
             particles_get_num_part_all(this%particles,stat_info_sub)
      
        
        !----------------------------------------------------
        ! Get all particles' (including ghosts) positions,
	! to build neighbor list.
	!----------------------------------------------------
        
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif

        CALL particles_get_x(this%particles,t_x, &
             num_part_all,stat_info_sub)
        CALL technique_build_list(this%tech,t_x, &
             num_part_all, symmetry,stat_info_sub)
        
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif   

        IF ( stat_info_sub /= 0 ) THEN
           PRINT *,"marching_integrate_VV : ", &
                "Building list failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
	! Compute mass density/number density for particles.
	!----------------------------------------------------

#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif
        
        CALL particles_compute_density(this%particles, &
             stat_info_sub)

#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif   
        
        IF ( stat_info_sub /= 0 ) THEN 
           PRINT *,"marching_integrate_VV : ", & 
                "Computing density failed!"
           stat_info = -1 
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! For symmtery calculation :
        ! swap the send and receive buffer,
        ! in order to send contribution of density
        ! of ghost particles to their host processes and
        ! receive contribution from other processes.
        !
        ! Even using single process, this has to be done,
        ! since it guarantees the boundary condition also.
        !----------------------------------------------------

#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif
        
        IF ( symmetry ) THEN
           
           CALL particles_map_ghost_put(this%particles, &
                l_map_x   = .TRUE., l_map_rho = .TRUE., &
                stat_info = stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_integrate_VV: ", &
                   "Receiving density from ghosts failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif   
        
        !-------------------------------------------------
        ! After computing density, the list of colloid
        ! boundary particles is created.
        ! Set colloidal boundary particles velocity
        ! accordingt to its translation and rotation speed.
        !-------------------------------------------------

        IF ( num_colloid > 0 ) THEN
           
           CALL particles_set_colloid_velocity(this%particles,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'marching_integrate_VV: ',&
                   'Setting colloid velocity failed !'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF

        !-------------------------------------------------
        ! After computing density, the list of wall
        ! boundary particles is created.
        ! Set wall boundary particles velocity
        ! accordingt to its translation speed.
        !-------------------------------------------------

        IF ( num_shear > 0 ) THEN

           CALL particles_set_boundary_velocity(this%particles,&
                stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_integrate_VV: ", &
                   "particles setting boundary failed !"
              stat_info = -1           
              GOTO 9999           
           END IF
           
        END IF

        !----------------------------------------------------
        ! After the density get contribution from ghosts
        ! on other proceses, create new ghost particles for 
        ! neigboring processes.
        !
        ! At this point, update ghosts with velocity
        ! which is needed for vgt and force.
        ! And ct which is needed for force.
        !
        ! Even using single process, this has to be done,
        ! since it guarantee the boundary condition also.
        !
	! Has to be done BEFORE building neighbor lists.
      	!----------------------------------------------------

#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif
        
        CALL particles_map_ghost_get(this%particles, &
             l_map_x   = .TRUE., l_map_rho = .TRUE., &
             l_map_v   = .TRUE., &
             l_map_ct  = (.NOT. Newtonian), &
             stat_info = stat_info_sub)
        
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif   

        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "marching_integrate_VV: ", &
                "Mapping ghosts with rho, v (,ct) failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! The number of real, all or ghost particles 
        ! should not be changed, shown here for clarity.        
        !----------------------------------------------------
        
        num_part_real  = &
             particles_get_num_part_real(this%particles,stat_info_sub)
        num_part_ghost = &
             particles_get_num_part_ghost(this%particles,stat_info_sub)
        num_part_all  = &
             particles_get_num_part_all(this%particles,stat_info_sub)
      
   
        !----------------------------------------------------
        ! After mapping velocity, according to different 
        ! boundary conditions, we have set up the veloicity 
        ! of boundary particles which are ghosts also.
        !----------------------------------------------------
        
        CALL particles_set_boundary_ghost_velocity(this%particles, &
             stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "marching_integrate_VV: ", &
                "Setting boundary ghosts velocity failed !"
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
           
           IF(stat_info_sub /=0) THEN
              PRINT *, "marching_integrate_VV: ", &
                   "Finding density extrem failed !"
              stat_info = -1
              GOTO 9999
           END IF

           CALL particles_set_stateEquation_rho_ref(this%particles,&
                stat_info_sub)
           
           IF(stat_info_sub /=0) THEN
              PRINT *, "marching_integrate_VV: ", &
                   "Setting density extreme failed !"
              stat_info = -1
              GOTO 9999
           END IF

        END IF

        
        !----------------------------------------------------
      	! Compute pressure for all particles, since ghosts 
        ! pressure is also needed to calculated force.
        !----------------------------------------------------

#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif
        
        CALL particles_compute_pressure(this%particles,&
             num_part_all,stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN
           PRINT *, "marching_integrate_VV: ", &
                "Computing pressure failed !" 
           stat_info = -1
           GOTO 9999
        END IF
        
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif   

        !----------------------------------------------------
      	! Compute pressure tensor for non-Newtonian case.
        ! (including ghost particles)
        !----------------------------------------------------
        
        IF  ( .NOT. Newtonian ) THEN
           
           CALL particles_compute_pressure_tensor(this%particles,&
                num_part_all,stat_info_sub)
           
           IF(stat_info_sub /=0) THEN
              PRINT *, "marching_integrate_VV: ", &
                   "Computing pressure tensor failed !"
              stat_info = -1
              GOTO 9999
           END IF

        END IF
        
        !----------------------------------------------------
        ! Compute interaction for all particles.
        ! (e.g. velocity gradient tensor
        ! (Non-Newotnian),
        ! pressure, viscous, thermal forces)
        !----------------------------------------------------

#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif
        
        CALL particles_compute_interaction(this%particles, &
             stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *,"marching_integrate_VV: ",&
                "Computing interaction failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif   
        
        !----------------------------------------------------
        ! For symmtery inter-process communication:
        !----------------------------------------------------

#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif
        
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
                l_map_f   = .TRUE., &
                l_map_fp  = .TRUE., &
                l_map_fv  = .TRUE., &
                l_map_fr  = Brownian, &
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
              PRINT *, "marching_integrate_VV: ", &
                   "Receiving force (stress, vgt, au) from ghosts failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! symmetry
        
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif   
        
        !----------------------------------------------------
	! Apply external / body force for real particles.
        !----------------------------------------------------
        
        CALL particles_apply_body_force(this%particles,&
             num_part_real,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *,"marching_integrate_VV: ", &
                "Applying body force failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! For non-Newtonian viscoelastic Oldroyd-B model.
        !----------------------------------------------------
        
        IF( .NOT. Newtonian ) THEN
           
           !-------------------------------------------------
           ! In case of egen-dynamics, compute accleration 
           ! of egenvalues and egenvectors, integrate eval
           ! and evec. Finally compute conformation tensor.
           !-------------------------------------------------
           
           IF ( eigen_dynamics ) THEN
              
              !----------------------------------------------
              ! Compute matrix element of eigen-dynamics 
              ! from velocity gradient tensor.
              !----------------------------------------------
            
              CALL particles_compute_evgt(this%particles, &
                   num_part_real,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV: ", &
                      "Computing evgt failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Compute accelerations of eigenvalues.
              !----------------------------------------------
          
              CALL particles_compute_aeval(this%particles,&
                   num_part_real,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV: ",&
                      "Computing aeval failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Compute accelerations of eigenvectors.
              !----------------------------------------------
            
              CALL particles_compute_aevec(this%particles, &
                   num_part_real,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV: ",&
                      "Computing aevec failed !"
                 stat_info = -1
                 GOTO 9999
              END IF

              !----------------------------------------------
              ! Integrate eigenvalues.
              !----------------------------------------------
              
              CALL particles_integrate_eval(this%particles, &
                   num_part_real,dt,0.5_MK,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *,&
                      "marching_integrate_VV: ", &
                      "Integrating eva 2nd time failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! Integrate eigenvectors.
              !----------------------------------------------
              
              CALL particles_integrate_evec(this%particles,&
                   num_part_real,dt,0.5_MK,stat_info_sub)
              
              IF (stat_info_sub /= 0 ) THEN
                 PRINT *,&
                      "marching_integrate_VV  ",&
                      "Integrating evec 2nd time failed !"
                 stat_info = -1
                 GOTO 9999
              END IF

              !----------------------------------------------
              ! Compute conformation tensor out of
              ! eigenvalues and eigenvectors second time.
              !
              ! In fact, for interaction(force) this is not 
              ! nessary, since eval and evec will be be
              ! integrated again, which will be used for
              ! computing ct before interaction(force).
              ! However, this can be used for output.
              !----------------------------------------------
              
              CALL particles_compute_ct(this%particles, &
                   num_part_real,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV : ", &
                      "Computing ct 2nd time failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              
           ELSE ! evolution of conformation tensor.
              
              CALL particles_compute_act(this%particles,&
                   num_part_real,stat_info_sub)
              
              IF (stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV : ",&
                      "Computing act failed !"
                 stat_info = -1
                 GOTO 9999              
              END IF

              !----------------------------------------------
              ! Integrate conformation tensor, using same 
              ! order as velocity second time.
              !----------------------------------------------
              
              CALL particles_integrate_ct(this%particles, &
                   num_part_real,dt,0.5_MK,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV : ", &
                      "Integrating ct 2nd time failed  !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF ! eigen-dynamics
           
        END IF ! non-Newtonian
        
#if 0        
        !----------------------------------------------------
        ! Since paire-wise interaction/force have been
        ! calculated, sum up the force on colloids,
        ! if colloids are present.
        !----------------------------------------------------

        
        IF ( num_colloid > 0 ) THEN
           
           !-------------------------------------------------
           ! Sum up force/torque exerted on parts of
           ! colloids on this process.
           !-------------------------------------------------
           
           CALL particles_collect_colloid_interaction(this%particles, &
                coll_drag,coll_torque,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_integrate_VV: ", &
                   "Summing up interaction on colloid locally has problem !"
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
              PRINT *, "marching_integrate_VV: ",&
                   "Summing up interaction on colloid globally has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Sum up force/torque exerted on colloids
           ! from all processes.
           !-------------------------------------------------
      
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time0 = &
                debug_get_time(global_debug,stat_info_sub)
        END IF
#endif     
           CALL colloid_compute_interaction(colloids,comm, &
                MPI_PREC,wall_drag_c(1:num_dim,1:num_dim*2),&
                stat_info_sub)
           
#ifdef __DEBUG_INTEGRATE_VV
        IF ( debug_flag == 3 ) THEN
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_index = debug_index + 1
           debug_time_record(debug_index) = &
                debug_time1 - debug_time0
        END IF
#endif   
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_integrate_VV: ", &
                   "Summing up interaction on colloid globally has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Apply body force on colloids.
           !-------------------------------------------------
        
           CALL colloid_apply_body_force(colloids,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_integrate_VV: ", &
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
              PRINT *, "marching_integrate_VV: ",&
                   "Computing colloids accelerations has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF !  num_colloid > 0
#endif
        
        !----------------------------------------------------
        ! If there is wall using symmetry, solid boundary
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
              PRINT *, "marching_integrate_VV: ", &
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
                MPI_PREC,wall_drag_p(1:num_dim,1:num_dim*2), &
                stat_info=stat_info_sub)
#endif     
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_integrate_VV: ",&
                   "Summing up particles contribution on boundary has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           CAll boundary_collect_colloid_interaction(tboundary,comm, &
                MPI_PREC, wall_drag_c(1:num_dim,1:num_dim*2),stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_integrate_VV: ",&
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
              PRINT *, "marching_integrate_VV: ",&
                   "Resetting boundary particles interaction failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
           CALL particles_reset_boundary_ghost_interaction(this%particles, &
                stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_integate_VV: ", &
                   "Resetting boundary interaction failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
           
        END IF ! num_shear > 0
        
        !----------------------------------------------------
	! Calculate new velocity 
        ! v(t+dt) =  v(t) + 0.5*(f(t)+f(t+dt)) i.e.,
        ! v(t+dt) =  v'(t+dt) + 0.5* f(t+dt) * dt 
	!----------------------------------------------------
        
        CALL particles_integrate_velocity(this%particles, &
             num_part_real,dt,0.5_MK,stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN
           PRINT *, "marching_integrate_VV: ", &
                "Updating velocity second time has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
#if 0        
        !----------------------------------------------------
        ! Adjust real particles' r/v after motion,
        ! according to boundary conditions,
	! in case they go out of the physical domain.
	!----------------------------------------------------
        CALL particles_adjust_particles(this%particles,&
             num_part_real,stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "marching_integrate_VV: ", &
                "Adjusting r or v failed ! "
           stat_info = -1
           GOTO 9999
        END IF      
#endif    
        
        !----------------------------------------------------
        ! Check if potential energy is needed.
        !----------------------------------------------------
        
        IF ( p_energy ) THEN
           
           CALL particles_integrate_potential_energy(this%particles, &
                num_part_real,dt,0.5_MK,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "marching_integrate_VV: ",&
                   "Updating potential energy second time has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! p_energy
    
#if 0

        !----------------------------------------------------
        ! If there are colloids, update their velocities
        ! as rigid bodies, using same integrator.
	!----------------------------------------------------
        
        IF ( num_colloid > 0 ) THEN
           
           CALL colloid_integrate_velocity(colloids,dt, &
                0.5_MK,stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_integrate_VV : ",&
                   "Updating colloids velocity second tinme failed !"
              stat_info = -1
              GOTO 9999
           END IF
        
           !-------------------------------------------------
           ! Set colloidal boundary particles velocity
           ! accordingt to its translation and rotation speed.
           !-------------------------------------------------
           
           CALL particles_set_colloid_velocity(this%particles,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'marching_integrate_VV: ',&
                   'Setting colloid velocity failed !'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! num_colloid > 0
        
#endif
        !----------------------------------------------------
        ! Update boundary :
        !
        ! such as walls' velocity in case of oscillating 
        ! shear;
        ! Or sheared length, in case of Lees-Edwards boundary.
        !----------------------------------------------------
        
        IF ( num_shear > 0 ) THEN
           
           CALL boundary_update_boundary(tboundary,&
                time+dt,time+dt,stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_integrate_VV: ", &
                   "Updating boundary failed  !"
              stat_info = -1
              GOTO 9999
           END IF
           
           CALL particles_set_boundary_velocity(this%particles,&
                stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, "marching_integrate_VV: ", &
                   "particles setting boundary failed !"
              stat_info = -1           
              GOTO 9999           
           END IF
           
        END IF ! num_shear > 0
        
#ifdef __DEBUG_INTEGRATE_VV
        
        SELECT CASE (debug_flag)
           
        CASE (2)
           CALL debug_substop(global_debug,rank,&
                "marching_integrate_VV",&
                debug_time_start,stat_info_sub)
        CASE (3)
           debug_time1 = &
                debug_get_time(global_debug,stat_info_sub)
           debug_time_record(1) = &
                debug_time1 - debug_time_start
           
           IF ( rank == 0 ) THEN
              
              CALL debug_write_time(global_debug,&
                   rank,debug_time_record(1:debug_index),stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "marching_integrate_VV: ", &
                      "debug_write_time has problem ! "
                 stat_info_sub = -1
                 GOTO 9999
              END IF
              
           END IF
           
        END SELECT
        
#endif 
        
9999	CONTINUE
        
        !----------------------------------------------------
        ! Release dynamic memories.
	!----------------------------------------------------
        
        IF (ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF
    
        IF (ASSOCIATED(t_x)) THEN
           DEALLOCATE(t_x)           
        END IF
        
        RETURN
        
      END SUBROUTINE marching_integrate_VV
      
      
      
