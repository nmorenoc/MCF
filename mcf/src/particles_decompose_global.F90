      SUBROUTINE particles_decompose_global(this, &
           l_map_x,    l_map_v,    &
           l_map_rho,  l_map_m,    &
           l_map_id,   l_map_f,    &
           l_map_vgt,  l_map_evgt, &
           l_map_eval, l_map_aeval,&
           l_map_evec, l_map_aevec,&
           l_map_ct,   l_map_act,  &
           l_map_u,    l_map_au, stat_info)
        !----------------------------------------------------
        !  Program      :   particles_decompose_global
        !----------------------------------------------------
        !
        !  Purpose      :   Broadcast variables from root
        !                   process to all processes.
        !                   What to broadcast depends on
        !                   the arguments.
        !
        !  Reference    :
        !
        !  Remark       : Positions are always needed for 
        !                 decomposing, therefore it doesn't 
        !                 need to specify anymore.
        !                
        !
        !  Revisions    :   V0.1 15.07.2009, original version.
        !
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        !       Dr. Marco Ellero's Emmy Noether Group,
        !       Prof. Dr. N. Adams' Chair of Aerodynamics,
        !	Faculty of Mechanical Engineering,
        !	Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Modules :
        !----------------------------------------------------
        
        USE ppm_module_map
        USE ppm_module_topo 
        USE ppm_module_find_duplicates
        
        !----------------------------------------------------
        ! Arguments.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)  :: this
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_x        
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_v
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_rho
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_m
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_id
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_f
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_vgt
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_evgt
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_eval
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_aeval
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_evec
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_aevec
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_ct
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_act
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_u
        LOGICAL, INTENT(IN), OPTIONAL   :: l_map_au     
        INTEGER, INTENT(OUT)            :: stat_info
        
        !--------------------------------
        ! Local variables
        !--------------------------------
        
        INTEGER                         :: stat_info_sub
	
        LOGICAL                         ::  map_x        
        LOGICAL                         ::  map_v
        LOGICAL                         ::  map_rho
        LOGICAL                         ::  map_m
        LOGICAL                         ::  map_id
        LOGICAL                         ::  map_f
        LOGICAL                         ::  map_vgt
        LOGICAL                         ::  map_evgt
        LOGICAL                         ::  map_eval
        LOGICAL                         ::  map_aeval
        LOGICAL                         ::  map_evec
        LOGICAL                         ::  map_aevec
        LOGICAL                         ::  map_ct
        LOGICAL                         ::  map_act
        LOGICAL                         ::  map_u
        LOGICAL                         ::  map_au

        INTEGER                         :: num_dim
        
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: num_wall_solid
        INTEGER                         :: num_part_wall_solid

        TYPE(Colloid), POINTER          :: colloids
        INTEGER                         :: num_colloid
        INTEGER, DIMENSION(:), POINTER  :: coll_num_physical_part
        INTEGER, DIMENSION(:), POINTER  :: coll_num_numerical_part
        
        INTEGER                         :: topo_id
        INTEGER                         :: map_type
        LOGICAL                         :: lcheck
        INTEGER				:: nid
        INTEGER, DIMENSION(:,:), POINTER:: ide
        INTEGER                         :: rank,comm,MPI_PREC

       
        !--------------------------------
        ! Initialization of variables.
        !--------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        map_x     = .FALSE.
        map_v     = .FALSE.
        map_rho   = .FALSE.
        map_m     = .FALSE.
        map_id    = .FALSE.
        map_f     = .FALSE.
        map_vgt   = .FALSE.
        map_evgt  = .FALSE.
        map_eval  = .FALSE.
        map_aeval = .FALSE.
        map_evec  = .FALSE.
        map_aevec = .FALSE.
        map_ct  = .FALSE.
        map_act = .FALSE.
        map_u   = .FALSE.
        map_au  = .FALSE.
        
        NULLIFY(tboundary)
        NULLIFY(colloids)
        NULLIFY(coll_num_physical_part)
        NULLIFY(coll_num_numerical_part)
        
        lcheck  = .TRUE.
        nid     = 0
      
        !--------------------------------
        ! Decide what to map.
        !--------------------------------
        
        IF( PRESENT(l_map_x) ) THEN	   
           map_x = l_map_x           
        END IF
        IF( PRESENT(l_map_v) ) THEN	   
           map_v = l_map_v           
        END IF
        IF( PRESENT(l_map_rho) ) THEN	   
           map_rho = l_map_rho           
        END IF
        IF( PRESENT(l_map_m) ) THEN	   
           map_m = l_map_m           
        END IF
        IF( PRESENT(l_map_id) ) THEN	   
           map_id = l_map_id           
        END IF
        IF( PRESENT(l_map_f) ) THEN	   
           map_f = l_map_f
        END IF
        IF( PRESENT(l_map_vgt) ) THEN
           map_vgt = l_map_vgt           
        END IF
        IF( PRESENT(l_map_evgt) ) THEN
           map_evgt = l_map_evgt           
        END IF
        IF( PRESENT(l_map_eval) ) THEN
           map_eval = l_map_eval           
        END IF
        IF( PRESENT(l_map_aeval) ) THEN
           map_aeval = l_map_aeval           
        END IF
        IF( PRESENT(l_map_evec) ) THEN
           map_evec = l_map_evec           
        END IF
        IF( PRESENT(l_map_aevec) ) THEN
           map_aevec = l_map_aevec           
        END IF
        IF( PRESENT(l_map_ct) ) THEN
           map_ct = l_map_ct           
        END IF
        IF( PRESENT(l_map_act) ) THEN
           map_act = l_map_act           
        END IF
        IF( PRESENT(l_map_u) ) THEN
           map_u = l_map_u           
        END IF
        IF( PRESENT(l_map_au) ) THEN
           map_au = l_map_au           
        END IF
        
        num_dim = this%num_dim
        
        rank    = &
             technique_get_rank(this%tech,stat_info_sub)
        topo_id = &
             technique_get_topo_id(this%tech,stat_info_sub)
        
     
#ifdef __MPI
        
        comm     = &
             technique_get_comm(this%tech,stat_info_sub)
        MPI_PREC = &
             technique_get_MPI_PREC(this%tech,stat_info_sub)
        
        
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        
        num_wall_solid = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        
        IF ( num_wall_solid > 0 ) THEN
           
           num_part_wall_solid = &
                boundary_get_num_part_wall_solid(tboundary,stat_info_sub)
           
           CALL MPI_BCAST(num_part_wall_solid, &
                1,MPI_INTEGER,&
                0,comm,stat_info_sub)  
           
           CALL boundary_set_num_part_wall_solid(tboundary,&
                num_part_wall_solid,stat_info_sub)
           
        END IF ! num_wall_solid > 0
        
        num_colloid  = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF( num_colloid >0 ) THEN
           
           !-------------------------------------------------
           ! If there are colloids,
           ! since the other parameters related to colloid
           ! have been read from file in parallel already,
           ! the num_physical_part, num_numerical_part are
           ! known only after we generate particles.
           ! Needed to be broadcasted to all the processes.
           !-------------------------------------------------
           
           CALL physics_get_colloid(this%phys,colloids,stat_info_sub)
           
           CALL colloid_get_num_physical_part(colloids,&
                coll_num_physical_part,stat_info_sub)        
           CALL colloid_get_num_numerical_part(colloids,&
                coll_num_numerical_part,stat_info_sub)


           CALL MPI_BCAST(coll_num_physical_part,&
                num_colloid,MPI_INTEGER,&
                0,comm,stat_info_sub)        
           CALL MPI_BCAST(coll_num_numerical_part,&
                num_colloid,MPI_INTEGER,&
                0,comm,stat_info_sub)
           
           
           CALL colloid_set_num_physical_part(colloids,&
                coll_num_physical_part,stat_info_sub)
           CALL colloid_set_num_numerical_part(colloids,&
                coll_num_numerical_part,stat_info_sub) 
           
        END IF
           
#endif
        
	!------------------------------------------
      	!  Map the particles onto the topology 
      	!------------------------------------------
        
        map_type = ppm_param_map_global
        
        !------------
        ! Position
        !------------
        
        IF ( map_x ) THEN
           
           CALL ppm_map_part(this%x, num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)        
           
           IF (stat_info_sub /= 0) THEN           
              PRINT *, "particles_decomposition_global : ",&
                   'Failed to start global mapping.'           
              stat_info = -1
              GOTO 9999           
           END IF
           
        END IF
        
        !----------------------
        ! Push to buffer
        !----------------------
        
        map_type = ppm_param_map_push
        
        
        !------------
        ! Velocity
        !------------        
        
        IF( map_v) THEN
           
           CALL ppm_map_part(this%v,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
	
        !------------        
        ! Density
        !------------        
        
        IF( map_rho) THEN
           
           CALL ppm_map_part(this%rho,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub) 
        END IF
        
        !------------        
        ! Mass
        !------------
        
        IF( map_m) THEN
           
           CALL ppm_map_part(this%m,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
        END IF
        
        !------------        
        ! IDs
        !------------
        
        IF (map_id) THEN
           
           CALL ppm_map_part(this%id,this%num_id,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
        END IF
       
        !------------
        ! Force
        !------------        
        
        IF( map_f) THEN
           
           CALL ppm_map_part(this%f,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF

        !--------------------------------
        ! Velocity Gradient tensor.
        !--------------------------------
        
        IF( map_vgt) THEN
           
           CALL ppm_map_part(this%vgt,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF

        !--------------------------------
        ! Egenvaector dynamics 
        ! Velocity Gradient tensor.
        !--------------------------------
        
        IF( map_evgt) THEN
           
           CALL ppm_map_part(this%evgt,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !--------------------------------
        ! Egenvaector dynamics 
        ! egenvalues.
        !--------------------------------
        
        IF( map_eval) THEN
           
           CALL ppm_map_part(this%eval,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF

        !--------------------------------
        ! Egenvaector dynamics 
        ! accleration of egenvalues.
        !--------------------------------
        
        IF( map_aeval) THEN
           
           CALL ppm_map_part(this%aeval,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !--------------------------------
        ! Egenvaector dynamics 
        ! egenvectors.
        !--------------------------------
        
        IF( map_evec) THEN
           
           CALL ppm_map_part(this%evec,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF

        !--------------------------------
        ! Egenvaector dynamics 
        ! accleration of egenvectors.
        !--------------------------------
        
        IF( map_aevec) THEN
           
           CALL ppm_map_part(this%aevec,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF

        !----------------------
        ! Conformation tensor.
        !----------------------
        
        IF( map_ct) THEN
           
           CALL ppm_map_part(this%ct,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !------------------------------------------
        ! Acceleration of Conformation tensor.
        !------------------------------------------
        
        IF( map_act) THEN
           
           CALL ppm_map_part(this%act,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !----------------------
        ! potential energy
        !----------------------
        
        IF( map_u) THEN
           
           CALL ppm_map_part(this%u,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
        END IF
        
        !--------------------------------
        ! potential energy accleration
        !--------------------------------
        
        IF( map_au) THEN
           
           CALL ppm_map_part(this%au,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
        END IF
        
        IF (stat_info_sub /= 0) THEN           
           PRINT *, &
                'particles_decomposition_global : ',&
                'Failed to push!'           
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !----------------------
        ! Send the data
        !----------------------        
        
        map_type = ppm_param_map_send
        
        CALL ppm_map_part(this%x,num_dim,&
             this%num_part_real,this%num_part_all,&
             topo_id,map_type,stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN           
           PRINT *, 'particles_decomposition_global : ',&
                'Failed to send positions!'
           stat_info = -1
           GOTO 9999           
        END IF
        
        !------------
        ! Pop
        !------------
        
        map_type = ppm_param_map_pop
        
        !--------------------------------
        ! Potential energy accleration
        !--------------------------------
        
        IF( map_au) THEN
           
           CALL ppm_map_part(this%au,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF

        !----------------------
        ! Potential energy
        !----------------------
        
        IF( map_u) THEN
           
           CALL ppm_map_part(this%u,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !------------------------------------------
        ! Acceleration of Conformation tensor.
        !------------------------------------------
        
        IF( map_act) THEN
           
           CALL ppm_map_part(this%act,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !----------------------
        ! Conformation tensor.
        !----------------------
        
        IF( map_ct) THEN
           
           CALL ppm_map_part(this%ct,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !--------------------------------
        ! Egenvaector dynamics 
        ! accleration of egenvectors.
        !--------------------------------
        
        IF( map_aevec) THEN
           
           CALL ppm_map_part(this%aevec,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !--------------------------------
        ! Egenvaector dynamics 
        ! egenvectors.
        !--------------------------------
        
        IF( map_evec) THEN
           
           CALL ppm_map_part(this%evec,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF

        !--------------------------------
        ! Egenvaector dynamics 
        ! accleration of egenvalues.
        !--------------------------------
        
        IF( map_aeval) THEN
           
           CALL ppm_map_part(this%aeval,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !--------------------------------
        ! Egenvaector dynamics 
        ! egenvalues.
        !--------------------------------
        
        IF( map_eval) THEN
           
           CALL ppm_map_part(this%eval,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !--------------------------------
        ! Egenvaector dynamics 
        ! Velocity Gradient tensor.
        !--------------------------------
        
        IF( map_evgt) THEN
           
           CALL ppm_map_part(this%evgt,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF

        !--------------------------------
        ! Velocity Gradient tensor.
        !--------------------------------
        
        IF( map_vgt) THEN
           
           CALL ppm_map_part(this%vgt,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        
        !------------
        ! Force
        !------------        
        
        IF( map_f) THEN
           
           CALL ppm_map_part(this%f,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !------------
        ! IDs
        !------------
        
        IF( map_id) THEN
           
           CALL ppm_map_part(this%id,this%num_id,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !------------        
        ! Mass
        !------------
        
        IF( map_m) THEN
           
           CALL ppm_map_part(this%m,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !------------
        ! Density
        !------------
        
        IF( map_rho) THEN
           
           CALL ppm_map_part(this%rho,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !------------        
        ! Velocity
        !------------
        
        IF( map_v) THEN
           
           CALL ppm_map_part(this%v,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
        END IF
        
        !------------
        ! Position
        !------------
        
        CALL ppm_map_part(this%x,num_dim,&
             this%num_part_real,this%num_part_all,&
             topo_id,map_type,stat_info_sub)
        
        
        this%num_part_real       = this%num_part_all
        this%num_part_ghost      = 0
        this%num_part_sym        = 0
        this%num_part_wall_sym   = 0
        this%num_part_wall_solid = 0
        this%num_part_le    = 0
        this%num_part_colloid    = 0

        
        IF (stat_info_sub /= 0) THEN
           PRINT *, &
                'particles_decomposition_global : ',&
                'Failed to pop !'
           stat_info = -1
           GOTO 9999
        END IF
        
        
      	!---------------------------------------------------
      	!  Check if particles have been mapped correctly.
      	!----------------------------------------------------
        
        CALL ppm_topo_check(xp=this%x,&
             NPART=this%num_part_real,&
             topo_ok=lcheck, info=stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN           
           PRINT *, &
                'particles_decomposition_global : ',&
                'Failed to check topology.'           
           stat_info = -1
           GOTO 9999           
        END IF
        
        
        IF (.NOT.lcheck) THEN
           PRINT *,'particles_decomposition_global : ',&
                'Particles not mapped correctly!'
           stat_info = -1
           GOTO 9999           
        END IF
        
      	!------------------------------------------
      	! Check if all particles are unique.
      	!------------------------------------------
        
        CALL ppm_find_duplicates(this%x,&
             num_dim,this%num_part_real,&
             nid,ide,stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN
           PRINT *, &
                'particles_decomposition_global : ',&
                'Failed to check particles!'
           stat_info = -1
           GOTO 9999
        END IF
        
        IF (nid /= 0) THEN
           PRINT *,'particles_decomposition_global : ',&
                'Found collocating particles!'
           PRINT *, "nid : ", nid
           PRINT *, "ide : ", ide(:,:)
           stat_info = -1
           GOTO 9999
        END IF
	
        
9999    CONTINUE
        
        IF(ASSOCIATED(coll_num_physical_part)) THEN
           DEALLOCATE(coll_num_physical_part)
        END IF
        
        IF(ASSOCIATED(coll_num_numerical_part)) THEN
           DEALLOCATE(coll_num_numerical_part)
        END IF
        
        IF(ASSOCIATED(ide)) THEN	   
           DEALLOCATE(ide,STAT=stat_info)	   
        END IF
        
        RETURN
        
      END SUBROUTINE particles_decompose_global
