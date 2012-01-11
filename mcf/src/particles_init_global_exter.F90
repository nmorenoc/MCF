      SUBROUTINE particles_init_global_exter(this,&
           d_x,d_v,d_rho,d_m,d_id, d_num_part,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_init_particles_global_exter
        !----------------------------------------------------
        !
        ! Purpose     :  When the particles' configuration
        !                are given externally,
        !                such as, defined by preprocessing,
        !                or from restart file of last run,
        !                this routine is called to set the 
        !                initial particles' configuration.
        !
        ! Routines    :
        !
        ! Remarks     :
        !
        ! References  :
        !
        ! Revisions   : V0.2 04.12 2009, check the work
        !               flow and supply with more comments.
        !            
        !               V0.1 20.03 2009, original version.
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
        ! Modules
        !----------------------------------------------------
        
        USE ppm_module_find_duplicates
        
        !----------------------------------------------------
        ! Arguments
    	!----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)  :: this
        REAL(MK), DIMENSION(:,:)        :: d_x
        REAL(MK), DIMENSION(:,:)        :: d_v
        REAL(MK), DIMENSION(:)          :: d_rho
        REAL(MK), DIMENSION(:)          :: d_m
        INTEGER, DIMENSION(:,:)         :: d_id
        INTEGER, INTENT(IN)             :: d_num_part
        INTEGER, INTENT(OUT)	        :: stat_info
        
        !----------------------------------------------------
    	! Local variables
    	!----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: p_energy
        INTEGER                         :: num_species
        INTEGER                         :: num_dim
        INTEGER, DIMENSION(2)           :: dim_x
        INTEGER, DIMENSION(2)           :: dim_v
        INTEGER, DIMENSION(2)           :: dim_id         
        INTEGER                         :: dim_rho
        INTEGER                         :: dim_m        
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: num_colloid
        TYPE(COLLOID), POINTER          :: colloids
        INTEGER, DIMENSION(:), POINTER  :: coll_num_numerical_part
        INTEGER                         :: sid,i

        
        INTEGER,DIMENSION(:,:),POINTER  :: ide
        INTEGER                         :: nid
    
        
#ifdef __DEBUG
        INTEGER                         :: debug_flag
        INTEGER                         :: debug_threshold
        REAL(MK)                        :: time_routine_start
        INTEGER                         :: rank
#endif     
        
        !----------------------------------------------------
        ! Initialization of variables.
    	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(bcdef)
        NULLIFY(tboundary)
        
        NULLIFY(colloids)
        NULLIFY(coll_num_numerical_part)
        
        
        NULLIFY(ide)

        
#ifdef __DEBUG
        debug_flag      = &
             debug_get_flag(global_debug,stat_info_sub)
        debug_threshold = 2
        rank            = &
             technique_get_rank(this%tech,stat_info_sub)
        
        IF(debug_flag >1 .OR. &
             debug_flag > debug_threshold) THEN
           CALL debug_substart(global_debug,rank, &
                'particles_init_particles_exter', &
                time_routine_start,stat_info_sub)
        END IF
#endif
        
        !---------------------------------------------------
        ! Contol parameters :
        !---------------------------------------------------
        
        p_energy = &
             control_get_p_energy(this%ctrl,stat_info_sub)
        
        !---------------------------------------------------
        ! Physics parameters :
        !
        ! Get number of species, dimension and colloids.
        ! A object of Colloid Class.
        !----------------------------------------------------
        
        num_species = &
             physics_get_num_species(this%phys,stat_info_sub) 
        num_dim     = this%num_dim
        
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        
        num_colloid = physics_get_num_colloid(this%phys,stat_info_sub)
        CALL physics_get_colloid(this%phys,colloids,stat_info_sub)
        
        
        !----------------------------------------------------
        ! number of all particles is
        ! the same as number of real particles,
        ! which is d_num_part from external file. 
        !----------------------------------------------------
        
        this%num_part_real = d_num_part
        this%num_part_all  = d_num_part
        CALL physics_set_num_part_tot(this%phys,&
             d_num_part,stat_info_sub)

        
        !----------------------------------------------------
        ! Get the dimensions of all input data,
        ! in order to allocate memory for
        ! a particle object.
        !----------------------------------------------------

        dim_x(1)  = SIZE(d_x,1)
        dim_x(2)  = SIZE(d_x,2)
        dim_v(1)  = SIZE(d_v,1)
        dim_v(2)  = SIZE(d_v,2)
        dim_rho   = SIZE(d_rho,1)
        dim_m     = SIZE(d_m,1)
        dim_id(1) = SIZE(d_id,1)
        dim_id(2) = SIZE(d_id,2)

        !----------------------------------------------------
        ! Check if the dimensions of data matches.
        !----------------------------------------------------
        
        IF ( dim_x(1)  /= num_dim .OR. &
             dim_x(2)  /= d_num_part ) THEN
           PRINT *, "particles_init_particles_exter : ", &
                "position dimensions don't match !"
           stat_info = -1
           GOTO 9999
        END IF

        IF ( dim_v(1)  /= num_dim .OR. &
             dim_v(2)  /= d_num_part ) THEN
           PRINT *, "particles_init_particles_exter : ", &
                "velocity dimensions don't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( dim_rho  /= d_num_part ) THEN
           PRINT *, "particles_init_particles_exter : ", &
                "rho dimension doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( dim_m  /= d_num_part ) THEN
           PRINT *, "particles_init_particles_exter : ", &
                "mass dimension doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF

        IF ( dim_id(1)  /= this%num_id .OR. &
             dim_id(2)  /= d_num_part ) THEN
           PRINT *, "particles_init_particles_exter : ", &
                "IDs dimensions don't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
    	!  Allocate memory for particles :  
        !
	!  x   : position
	!  v   : velocity
        !  rho : density
        !  m   : mass
	!  IDs : particle id,  species id
    	!----------------------------------------------------
        
        IF (ASSOCIATED(this%x)) THEN
           DEALLOCATE(this%x,STAT=stat_info_sub)
        END IF
        
        IF (ASSOCIATED(this%v)) THEN
           DEALLOCATE(this%v,STAT=stat_info_sub)
        END IF
        
        IF (ASSOCIATED(this%rho)) THEN
           DEALLOCATE(this%rho,STAT=stat_info_sub)
        END IF
        
        IF (ASSOCIATED(this%m)) THEN
           DEALLOCATE(this%m,STAT=stat_info_sub)
        END IF
        
        IF (ASSOCIATED(this%id)) THEN
           DEALLOCATE(this%id,STAT=stat_info_sub)
        END IF
        
        ALLOCATE(this%x(dim_x(1),dim_x(2)), STAT=stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "particles_init_particles_exter : " , &
                "Allocating x has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        ALLOCATE(this%v(dim_v(1),dim_v(2)), STAT=stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "particles_init_particles_exter : " , &
                "Allocating v has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        ALLOCATE(this%rho(dim_rho), STAT=stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "particles_init_particles_exter : " , &
                "Allocating rho has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        ALLOCATE(this%m(dim_m), STAT=stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "particles_init_particles_exter : " , &
                "Allocating x has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        ALLOCATE(this%id(dim_id(1),dim_id(2)), STAT=stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "particles_init_particles_exter : " , &
                "Allocating id has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! u  :  potential energy, initially 1.0.
        ! It is not always needed.
        !----------------------------------------------------
        
        IF( p_energy ) THEN
           
           IF (ASSOCIATED(this%u)) THEN
              DEALLOCATE(this%u,STAT=stat_info_sub)
           END IF
           
           ALLOCATE(this%u(d_num_part), STAT=stat_info_sub)
           
           this%u(:) = 1.0_MK
           
        END IF
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "particles_init_particles_exter : " , &
                "Allocating memory for members has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Copying data to particles' members.
        !----------------------------------------------------
        
        this%x(:,:)  = d_x(:,:)
        this%v(:,:)  = d_v(:,:)
        this%rho(:)  = d_rho(:)
        this%m(:)    = d_m(:)
        this%id(:,:) = d_id(:,:)
        
        !----------------------------------------------------
    	! Check if each particle is unique.
    	!----------------------------------------------------
        
        nid = 0    
        CALL ppm_find_duplicates(this%x,num_dim,d_num_part,&
             nid,ide,stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN
           PRINT *, "particles_init_particles_exter : ", &
                "Error by checking duiplicates !"
           stat_info = -1
           GOTO 9999           
           
        END IF
        
        IF ( nid > 0 ) THEN
           PRINT *, &
                "particles_init_particles_exter : ", &
                "Found collocating particles !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Count & Record :
        !
        ! 1 How many fluid particles.
        ! 2 How many symmmetry boundary particles.
        ! 3 How many wall using symmmetry boundary particles.
        ! 4 How many wall solid boundary particles.
        ! 5 How many Lees-Edwards boundary particles.
        ! 6 How many numerical particles
        !   constitute each colloid.
        ! 7 How many of colloid boundary particles
        !   in total.
        !----------------------------------------------------
        
        this%num_part_fluid      = 0
        this%num_part_sym        = 0
        this%num_part_wall_sym   = 0
        this%num_part_wall_solid       = 0
        this%num_part_wall_solid_real  = 0
        this%num_part_wall_solid_ghost = 0
        this%num_part_le         = 0
        this%num_part_colloid = 0
        
        IF( num_colloid > 0 ) THEN
           
           ALLOCATE(coll_num_numerical_part(num_colloid))
           
           coll_num_numerical_part(:) = 0
           
        END IF
        
        !----------------------------------------------------
        ! Loop over all particles.
        !----------------------------------------------------
        
        Do i = 1, d_num_part
           
           !-------------------------------------------------
           ! Get species ID.           
           !-------------------------------------------------
           
           sid = d_id(this%sid_idx,i)
           
           IF( sid < -2*num_dim .OR. &
                sid > num_colloid) THEN
              
              PRINT *, "particles_init_particles_exter : ",&
                   "Species ID is wrong !"
              stat_info = -1
              GOTO 9999
              
           ELSE IF ( sid == 0 ) THEN
              
              this%num_part_fluid = &
                   this%num_part_fluid + 1
              
              
           ELSE IF ( sid < 0 ) THEN
              
              !----------------------------------------------
              ! Can only be solid wall boundary particles.
              !----------------------------------------------
              
              this%num_part_wall_solid = &
                   this%num_part_wall_solid + 1
              
              
           ELSE IF( sid > 0 ) THEN
              
              coll_num_numerical_part(sid) = &
                   coll_num_numerical_part(sid) + 1
              
              this%num_part_colloid = &
                   this%num_part_colloid + 1
              
           END IF ! sid
           
        END DO ! i = 1, d_num_part
        
        !----------------------------------------------------
        ! Save the number of numerical solid wall 
        ! boundary particles.
        !----------------------------------------------------
        
        CALL boundary_set_num_part_wall_solid(tboundary, &
             this%num_part_wall_solid,stat_info_sub)
        
        !----------------------------------------------------
        ! Save the number of numerical particles
        ! for each colloid.
        !----------------------------------------------------
        
        IF ( num_colloid > 0 ) THEN
           
           CALL colloid_set_num_numerical_part(colloids,&
                coll_num_numerical_part,stat_info_sub)
           
           IF( stat_info_sub > 0 ) THEN
              PRINT *, "particles_init_particles_exter : ",&
                   "Setting numerical particles is wrong !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
#if 0
        !----------------------
        !  For debug.
        !----------------------
        CALL debug_write_output(global_debug,rank,&
             "particles_init_particles_exter",&
             "x_real",0,this%x,1,this%num_part_real,stat_info_sub)
        
#endif 
        
9999    CONTINUE        
        
        IF (ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF
        
        IF(ASSOCIATED(coll_num_numerical_part)) THEN
           DEALLOCATE(coll_num_numerical_part)
        END IF
        
        IF ( ASSOCIATED(ide) ) THEN
           DEALLOCATE(ide)
        END IF
        
#ifdef __DEBUG        
        IF(debug_flag >1 .OR. debug_flag > debug_threshold) THEN
           CALL debug_substop(global_debug,rank,&
                'particles_init_particles_exter',&
                time_routine_start,stat_info_sub)
        END IF
#endif        
        
        RETURN
        
      END SUBROUTINE particles_init_global_exter
      
