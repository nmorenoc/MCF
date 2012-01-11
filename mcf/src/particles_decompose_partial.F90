      SUBROUTINE particles_decompose_partial(this,&
           l_map_x,    l_map_v,     &
           l_map_rho,  l_map_m,     &
           l_map_id,   l_map_f,     &
           l_map_vgt,  l_map_evgt,  &
           l_map_eval, l_map_aeval, &
           l_map_evec, l_map_aevec, &
           l_map_ct,   l_map_act,   &
           l_map_u,    l_map_au,stat_info)
        !----------------------------------------------------
        ! Program     :  particles_decompose_partial
        !----------------------------------------------------
        !
        ! Purpose     :  Communicate variables between
        !                neibouring processes.
        !                What to comminicate depends on
        !                the arguments.
        !
        ! Reference   :
        !
        ! Remark      :  Positions are always needed for 
        !                decomposing at the moment.
        !                
        !
        ! Revisions   :  V0.3 12.04.2011, check collocating
        !                particels, but not treat them as
        !                error, since it could happen in
        !                the presence of dense colloidal
        !                suspension.
        !
        !                V0.2 06.10.2010, try to output the
        !                particles' positions, when it failed
        !                to decompose.
        !
        !                V0.1 15.07.2009, original version.
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
        ! Modules used.
        !----------------------------------------------------
        
        USE ppm_module_map
        USE ppm_module_topo_check
        USE ppm_module_find_duplicates
        
      	!----------------------------------------------------
      	!  Arguments
      	!----------------------------------------------------
        
        TYPE(Particles),INTENT(INOUT)   :: this
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
        
      	!----------------------------------------------------
      	! Local variables
      	!----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_dim
        INTEGER                         :: topo_id
        INTEGER                         :: map_type         
        INTEGER				:: nid
        INTEGER, DIMENSION(:,:), POINTER:: ide
        LOGICAL	 		 	:: lcheck         
        
        
	!----------------------------------------------------
      	! Local/MAPPING
      	!----------------------------------------------------
	
        LOGICAL                         :: map_x
        LOGICAL                         :: map_v
        LOGICAL                         :: map_rho
        LOGICAL                         :: map_m
        LOGICAL                         :: map_id
        LOGICAL                         :: map_f
        LOGICAL                         :: map_vgt
        LOGICAL                         :: map_evgt
        LOGICAL                         :: map_eval
        LOGICAL                         :: map_aeval
        LOGICAL                         :: map_evec
        LOGICAL                         :: map_aevec
        LOGICAL                         :: map_ct
        LOGICAL                         :: map_act
        LOGICAL                         :: map_u
        LOGICAL                         :: map_au
        
        INTEGER                         :: rank

#ifdef __DEBUG

        REAL(MK),DIMENSION(:,:),POINTER :: t_x

        NULLIFY(t_x)
#endif
        
      	!----------------------------------------------------
      	! Initialize variables
      	!----------------------------------------------------
        
        stat_info       = 0
        stat_info_sub   = 0
        
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
        map_ct    = .FALSE.
        map_act   = .FALSE.
        map_u     = .FALSE.
        map_au    = .FALSE.
        
        lcheck  = .TRUE.
        nid     = 0
        
      	!----------------------------------------------------
        ! Change flag of each variable according to input
        ! arguments. The flags indicate whether or not
        ! the variables needs to be communicated.
      	!----------------------------------------------------
        
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
        
        num_dim = physics_get_num_dim(this%phys,stat_info_sub)
        topo_id = technique_get_topo_id(this%tech,stat_info_sub)
        rank    = technique_get_rank(this%tech,stat_info_sub)
        
	!----------------------------------------------------
      	! Decomposition starts
      	!----------------------------------------------------


	!----------------------------------------------------
      	! Choose mapping partial
      	!----------------------------------------------------
        
        map_type = ppm_param_map_partial
        
        IF ( map_x ) THEN
           
           CALL ppm_map_part(this%x,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Starting decomposion x failed!'	   
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
	!----------------------------------------------------
     	! Push the data into buffer
     	!----------------------------------------------------
        
        map_type = ppm_param_map_push
        
        IF ( map_v ) THEN
           
           CALL ppm_map_part(this%v,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing v failed!'	   
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        IF ( map_rho) THEN
           
           CALL ppm_map_part(this%rho,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing rho failed !'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        IF ( map_m) THEN
           
           CALL ppm_map_part(this%m,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing m failed !'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        IF ( map_id) THEN
           
           CALL ppm_map_part(this%id,this%num_id,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing id failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        IF ( map_f) THEN
           
           CALL ppm_map_part(this%f,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing f failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF

        IF ( map_vgt) THEN
           
           CALL ppm_map_part(this%vgt,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing vgt failed!'
              stat_info = -1
              GOTO 9999
           END IF
        END IF
        
        !----------------------------------------------------
        ! Egenvaector dynamics Velocity Gradient tensor.
        !----------------------------------------------------
        
        IF( map_evgt) THEN
           
           CALL ppm_map_part(this%evgt,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing evgt failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Egenvaector dynamics egenvalues.
        !----------------------------------------------------
        
        IF( map_eval) THEN
           
           CALL ppm_map_part(this%eval,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing eval failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Egenvaector dynamics accleration of egenvalues.
        !----------------------------------------------------
        
        IF( map_aeval) THEN
           
           CALL ppm_map_part(this%aeval,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)

           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing aeval failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Egenvaector dynamics egenvectors.
        !----------------------------------------------------
        
        IF( map_evec) THEN
           
           CALL ppm_map_part(this%evec,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing evec failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF

        !----------------------------------------------------
        ! Egenvaector dynamics accleration of egenvectors.
        !----------------------------------------------------
        
        IF( map_aevec) THEN
           
           CALL ppm_map_part(this%aevec,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)

           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing aevec failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF

        !----------------------------------------------------
        ! Conformation tensor.
        !----------------------------------------------------
        
        IF( map_ct) THEN
           
           CALL ppm_map_part(this%ct,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing ct failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Acceleration of Conformation tensor.
        !----------------------------------------------------
        
        IF( map_act) THEN
           
           CALL ppm_map_part(this%act,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing act failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
        IF ( map_u) THEN
           
           CALL ppm_map_part(this%u,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing u failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        IF ( map_au) THEN
           
           CALL ppm_map_part(this%au,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Pushing au failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
	!----------------------------------------------------
     	! Send the buffer
     	!----------------------------------------------------
        
        map_type = ppm_param_map_send
        
        CALL ppm_map_part(this%x,num_dim,&
             this%num_part_real,this%num_part_all,&
             topo_id,map_type,stat_info_sub)     
        
        
        IF (stat_info_sub /=0 ) THEN
           PRINT *, 'particles_decomposition_partial : ', &
                'Sending data failed!'
           stat_info = -1
           GOTO 9999
        END IF
        
	!----------------------------------------------------
     	! Pop data from the buffer
     	!----------------------------------------------------
        
        map_type = ppm_param_map_pop

        IF ( map_au) THEN
           
           CALL ppm_map_part(this%au,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)	
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping au failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
        IF ( map_u) THEN
           
           CALL ppm_map_part(this%u,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)	
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping u failed!'
              stat_info = -1
              GOTO 9999
           END IF

        END IF
       
        
        IF( map_act) THEN
           
           CALL ppm_map_part(this%act,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping act failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF

        
        IF( map_ct) THEN
           
           CALL ppm_map_part(this%ct,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping ct failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Egenvaector dynamics accleration of egenvectors.
        !----------------------------------------------------
        
        IF( map_aevec) THEN
           
           CALL ppm_map_part(this%aevec,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping aevec failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Egenvaector dynamics egenvectors.
        !----------------------------------------------------
        
        IF( map_evec) THEN
           
           CALL ppm_map_part(this%evec,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping evec failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF

        !----------------------------------------------------
        ! Egenvaector dynamics accleration of egenvalues.
        !----------------------------------------------------
        
        IF( map_aeval) THEN
           
           CALL ppm_map_part(this%aeval,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping aeval failed!'
              stat_info = -1
              GOTO 9999
           END IF

        END IF
        
        !----------------------------------------------------
        ! Egenvaector dynamics egenvalues.
        !----------------------------------------------------
        
        
        IF( map_eval) THEN
           
           CALL ppm_map_part(this%eval,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping eval failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Egenvaector dynamics Velocity Gradient tensor.
        !----------------------------------------------------
        
        IF( map_evgt) THEN
           
           CALL ppm_map_part(this%evgt,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)

            IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping evgt failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF

        
        IF( map_vgt) THEN
           
           CALL ppm_map_part(this%vgt,num_dim**2,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping vgt failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
        IF( map_f) THEN
           
           CALL ppm_map_part(this%f,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub) 
        
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping f failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
        IF ( map_id) THEN
           
           CALL ppm_map_part(this%id,this%num_id,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping id failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
        IF ( map_m) THEN
           
           CALL ppm_map_part(this%m,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping m failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
        IF ( map_rho) THEN
           
           CALL ppm_map_part(this%rho,&
                this%num_part_real,this%num_part_all,&
                topo_id,map_type,stat_info_sub)	
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping rho failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
                

        IF ( map_v) THEN
           
           CALL ppm_map_part(this%v,num_dim,&
                this%num_part_real,this%num_part_all,&
                topo_id, map_type, stat_info_sub) 
           
            IF (stat_info_sub /= 0) THEN
              PRINT *, 'particles_decomposition_partial : ',&
                   'Poping v failed!'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
        CALL ppm_map_part(this%x,num_dim,&
             this%num_part_real,this%num_part_all,&
             topo_id,map_type,stat_info_sub)
        
        
        IF (stat_info_sub /=0 ) THEN
           PRINT *, 'particles_decomposition_partial : ',&
                'Poping x failed!'
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Update particles number after mapping
        !----------------------------------------------------
        
        this%num_part_real     = this%num_part_all 
        this%num_part_ghost    = 0
        this%num_part_fluid    = 0
        this%num_part_sym      = 0
        this%num_part_wall_sym = 0
        this%num_part_wall_solid_real  = 0
        this%num_part_wall_solid_ghost = 0
        this%num_part_colloid = 0

        
        !----------------------------------------------------
        ! Check if particles have been mapped correctly
        !----------------------------------------------------
          
        CALL ppm_topo_check(xp=this%x,&
             NPART=this%num_part_real,&
             topo_ok=lcheck, info=stat_info_sub)
        
        IF (stat_info_sub /= 0 ) THEN
           PRINT *, 'particles_decomposition_partial : ',&
                'ppm_topo_check failed!'
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( .NOT. lcheck ) THEN
           PRINT *,'particles_decomposition_partial : ',&
                'Particles not mapped correctly!'
           stat_info = -1
           GOTO 9999
        END IF
    
        
	!----------------------------------------------------
      	! Check if all particles are unique
      	!----------------------------------------------------
        
        CALL ppm_find_duplicates(this%x,num_dim,&
             this%num_part_real,&
             nid,ide,stat_info_sub)
	
        IF (stat_info_sub /= 0) THEN
           PRINT *,'particles_decomposition_partial : ',&
                'Failed to check particles !'
           stat_info = -1 
           GOTO 9999
        END IF
        
        IF (nid /= 0) THEN	   
           PRINT *,'particles_decomposition_partial : ',&
                'Found collocating particles !'
           PRINT *, "nid : ", nid
           PRINT *, "ide : ", ide(:,:)
           !stat_info = -1 
           !GOTO 9999
        END IF
        
        
9999    CONTINUE	
        
        IF( ASSOCIATED(ide) ) THEN	   
           DEALLOCATE(ide,STAT=stat_info_sub)	   
        END IF
        
        
#if __DEBUG
        
        !----------------------------------------------------
        ! For debugging, output x when it crashes.
        !----------------------------------------------------
        
        IF ( stat_info /= 0 ) THEN
           
           CALL particles_get_x(this,t_x, &
                this%num_part_real,stat_info_sub)
           CALL debug_write_output(global_debug,rank,"decompose", &
                'x_decom_real',0,t_x,1,this%num_part_real,stat_info_sub)
           
        END IF
        
        IF(ASSOCIATED(t_x)) THEN
           DEALLOCATE(t_x)
        END IF
        
#endif
        
        RETURN 
	
      END SUBROUTINE particles_decompose_partial
      
