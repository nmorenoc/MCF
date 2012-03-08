      SUBROUTINE io_read_physics_config(this,phys,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_read_physics_config
        !----------------------------------------------------
        !  Purpose    : Reading all physics parameters, 
        !               including colloid parameters,
        !               from physics config file.
        !
        !  Reference  : Ellis T.M.R. et al. 
        !               Fortran 90 Programming.
        !
        !  Remark     : 
        !
        !  Version    :  V0.1, 01.01.2009, original version.
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
        
        TYPE(IO), INTENT(IN)            :: this
        TYPE(Physics), INTENT(INOUT)    :: phys
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
      	! Local variables 
      	!----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                        	:: i,j,idx
        INTEGER                        	:: ilen,ios
        INTEGER                         :: iline,ilenphys
        CHARACTER(LEN=MAX_CHAR)        	:: cbuf,cvalue,carg
        LOGICAL                        	:: lExist
        
        !----------------------------------------------------
        ! Physics variables.
        !----------------------------------------------------
        
        INTEGER                         :: num_species
        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(3)          :: min_phys
        REAL(MK), DIMENSION(3)          :: max_phys
        INTEGER                         :: lattice_type
        INTEGER, DIMENSION(3)           :: num_part_dim
        REAL(MK)                        :: cut_off

        REAL(MK)                        :: dt
        INTEGER                         :: step_start
        INTEGER                         :: step_end        
        REAL(MK)                        :: time_start
        REAL(MK)                        :: time_end

        REAL(MK)                        :: rho
        REAL(MK)                        :: eta
        REAL(MK)                        :: eta_coef
        REAL(MK)                        :: ksai
        REAL(MK)                        :: kt
        REAL(MK)                        :: c
        REAL(MK)                        :: rho_ref
        REAL(MK)                        :: gamma
        
        INTEGER                         :: relax_type
        REAL(MK)                        :: dt_relax
        INTEGER                         :: step_relax
        REAL(MK)                        :: time_relax
        REAL(MK)                        :: disorder_level
        REAL(MK)                        :: kt_relax
        REAL(MK)                        :: c_relax

        REAL(MK)                        :: tau
	REAL(MK)                        :: tau_sm
	REAL(MK)                        :: k_sm
        REAL(MK)                        :: n_p
        REAL(MK)                        :: kt_p
        LOGICAL                         :: eigen_dynamics
        REAL(MK), DIMENSION(6)          :: eval
        REAL(MK)                        :: eval_tolerance
        REAL(MK), DIMENSION(6,6)        :: evec
        LOGICAL                         :: evec_normalize
        REAL(MK)                        :: evec_tolerance
              
        INTEGER                         :: body_force_type
        REAL(MK), DIMENSION(3)          :: body_force
        REAL(MK), DIMENSION(3)          :: body_force_d
        INTEGER                         :: flow_direction
        REAL(MK)                        :: flow_width
        REAL(MK)                        :: flow_v
        INTEGER                         :: flow_adjust_freq
         
        !----------------------------------------------------
        ! colloid variables.
        !----------------------------------------------------
        
        INTEGER                                 :: num_colloid
        REAL(MK)                                :: coll_rho
        INTEGER                                 :: coll_rho_type
        LOGICAL                                 :: coll_translate
        LOGICAL                                 :: coll_rotate
        INTEGER                                 :: coll_place
        INTEGER                                 :: coll_noslip
        INTEGER                                 :: coll_body_force_type
        REAL(MK), DIMENSION(3)                  :: coll_body_force
        
       
        REAL(MK)                                :: cc_lub_cut_off
        REAL(MK)                                :: cc_lub_cut_on
        REAL(MK)                                :: cc_repul_cut_off
        REAL(MK)                                :: cc_repul_cut_on
        REAL(MK)                                :: cc_repul_F0
        
        REAL(MK)                                :: cw_lub_cut_off
        REAL(MK)                                :: cw_lub_cut_on
        REAL(MK)                                :: cw_repul_cut_off
        REAL(MK)                                :: cw_repul_cut_on
        REAL(MK)                                :: cw_repul_F0


        TYPE(Colloid), POINTER                  :: colloids
        INTEGER,ALLOCATABLE,DIMENSION(:)        :: coll_shape
        REAL(MK),ALLOCATABLE, DIMENSION(:,:)    :: coll_radius
        INTEGER, ALLOCATABLE, DIMENSION(:)      :: coll_freq
        REAL(MK),ALLOCATABLE, DIMENSION(:)      :: coll_m
        REAL(MK),ALLOCATABLE, DIMENSION(:,:)    :: coll_mmi
        REAL(MK),ALLOCATABLE,DIMENSION(:,:)     :: coll_x        
        REAL(MK),ALLOCATABLE,DIMENSION(:,:)     :: coll_v
        REAL(MK),ALLOCATABLE, DIMENSION(:,:)    :: coll_acc_vector
        REAL(MK),ALLOCATABLE, DIMENSION(:,:)    :: coll_theta    
        REAL(MK),ALLOCATABLE,DIMENSION(:,:)     :: coll_omega
        INTEGER                                 :: coll_index
        INTEGER                                 :: shape_index
        INTEGER                                 :: radius_index
        INTEGER                                 :: freq_index
        INTEGER                                 :: m_index
        INTEGER                                 :: mmi_index
        INTEGER                                 :: x_index
        INTEGER                                 :: v_index
        INTEGER                                 :: rot_v_index
        INTEGER                                 :: theta_index
        INTEGER                                 :: omega_index
        
        !----------------------------------------------------
        ! boundary variables.
        !----------------------------------------------------
     
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER, DIMENSION(6)           :: bcdef
        INTEGER, DIMENSION(6)           :: shear_type
        REAL(MK), DIMENSION(3,3)        :: shear_rate
        REAL(MK), DIMENSION(3,6)        :: shear_v
        REAL(MK), DIMENSION(6)          :: shear_freq
        INTEGER                         :: wall_rho_type
        INTEGER                         :: wall_noslip
        
#ifdef __DEBUG
        INTEGER                         :: debug_threshold
#endif
        
	!----------------------------------------------------
      	! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        !----------------------------------------------------
        ! boundary or
        ! colloids shouldn't be released afterwards.
        !----------------------------------------------------
        
        min_phys(:)     = 0.0_MK
        max_phys(:)     = 0.0_MK
        num_part_dim(:) = 0
        body_force(:)   = 0.0_MK
        body_force_d(:) = 0.0_MK
        eigen_dynamics  = .FALSE.
        eval(:)         = 0.0_MK
        eval_tolerance  = mcf_machine_zero
        evec(:,:)       = 0.0_MK
        evec_normalize  = .FALSE.
        evec_tolerance  = mcf_machine_zero
        
        num_colloid   = 0
      
        NULLIFY(colloids)
        ALLOCATE(colloids)        
        coll_index   = 0
        
        bcdef(:)        = 0
        shear_type(:)   = 1
        shear_rate(:,:) = 0.0_MK
        shear_v(:,:)    = 0.0_MK
        shear_freq(:)   = 0.0_MK
        wall_rho_type   = 0
        wall_noslip     = 1
        NULLIFY(tboundary)
        ALLOCATE(tboundary)
        
        
#ifdef __DEBUG
        debug_threshold = 3
        IF ( debug_threshold < 3) THEN
           PRINT *, 'io_read_physics_config : starting'    
        END IF
#endif
        
        !----------------------------------------------------
        ! Check if the name config file empty.
      	!----------------------------------------------------
        
        ilenphys = LEN_TRIM(this%physics_config_file)
        
        IF (ilenphys < 1) THEN
           PRINT *,"io_read_physics_config : ",&
                "No physics config file name given !"
           stat_info = -1
           GOTO 9999
        END IF
        
      	!----------------------------------------------------
      	! Check if the physics config file exists.
      	!----------------------------------------------------
        
        INQUIRE(FILE=this%physics_config_file,EXIST=lExist)
	
        IF (.NOT.lExist) THEN
           WRITE(cbuf,'(2A)') &
                'No physics config file with name : ',&
                this%physics_config_file(1:ilenphys)
           PRINT *, 'io_read_physics_config : "', cbuf           
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Open physics config file for reading.
        !----------------------------------------------------
        
        OPEN(this%physics_config_unit,&
             FILE=this%physics_config_file,&
             IOSTAT=ios,ACTION='READ')
        
        IF (ios /= 0) THEN
           WRITE(cbuf,'(2A)')&
                'Failed to open file ',&
                this%physics_config_file(1:ilenphys)
           PRINT *,'io_read_physics_config : ', cbuf
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Scan the file line by line.
        !----------------------------------------------------
        
        iline = 0
	
      	DO 
          
          !--------------------------------------------------
          ! Increase line counter.
          !--------------------------------------------------
          
          iline = iline + 1
          
          !--------------------------------------------------
          ! Read one line.
          !--------------------------------------------------
          
          READ(this%physics_config_unit,'(A)',END=9999,ERR=200) cbuf 
          ilen = LEN_TRIM(cbuf)
          
#ifdef __DEBUG
          !--------------------------------------------------
          ! Extensive Debug: print each line of Ctrl file.
          !--------------------------------------------------
          IF ( debug_threshold < 2 ) THEN             
             PRINT *,'io_read_physics_config : ', cbuf(1:ilen)             
          END IF
#endif          
          
          !--------------------------------------------------
          ! Skip empty line or 
          ! comment line started with symbol #.
          !--------------------------------------------------
          
          IF (ilen < 1 .OR. cbuf(1:1) == '#' ) THEN
             CYCLE
          END IF
          
          !--------------------------------------------------
          ! Remove spaces of the line being read.
          !--------------------------------------------------
          
          j = 0
          DO i=1,ilen
             IF (cbuf(i:i) /=' ' .AND. cbuf(i:i) /= '\t' ) THEN
                j = j + 1
                cbuf(j:j) = cbuf(i:i)
             END IF
          END DO
          
          !--------------------------------------------------
          ! Update the length of this string.
          !--------------------------------------------------
          
          ilen = j
          
          !--------------------------------------------------
          ! After revoing the spaces,
          ! skip comment line started with symbol #.
          !--------------------------------------------------
          
          IF (cbuf(1:1) == '#' ) THEN
             CYCLE
          END IF
          
          !--------------------------------------------------
          ! Find the position of '='.
          !--------------------------------------------------
          
          idx = INDEX(cbuf,'=')
          
          !--------------------------------------------------
          ! Exit if '=' is missing !
          !--------------------------------------------------
          
          IF (idx < 0) THEN             
             WRITE(cbuf,'(A,I5)') &
                  "Incorrect line : ",&
                  iline, "No = found !"             
             Print *,'io_read_physics_config : ', cbuf             
             stat_info = -1
             GOTO 9999             
          END IF
          
          !--------------------------------------------------
          ! Get argument and its value.
          !--------------------------------------------------
          
          carg   = ADJUSTL(cbuf(1:idx-1))
          cvalue = ADJUSTL(cbuf(idx+1:ilen))
          
          !--------------------------------------------------
          ! Convert all arguments to capital letter.
          !--------------------------------------------------
          
          CALL tool_uppercase(this%io_tool,carg,idx-1,stat_info)
          
#ifdef __DEBUG
          IF ( debug_threshold < 1) THEN
             PRINT *, 'io_read_physics_config : ', carg
             PRINT *, 'io_read_physics_config : ', cvalue
          END IF
#endif   
          
          IF (carg == 'NUM_SPECIES') THEN
             
             !-----------------------------------------------
             ! Get number of species and
             ! set it into physics object.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) num_species
             
             CALL physics_set_num_species(phys,num_species,stat_info_sub)
             
             
          ELSE IF (carg == 'NUM_DIM') THEN
             
             !-----------------------------------------------
             ! Get space dimensionality and
             ! set it into physics object.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) num_dim
             
             CALL physics_set_num_dim(phys,num_dim,stat_info_sub)            
             
             
          ELSE  IF (carg == 'MIN_PHYS') THEN
             
             !-----------------------------------------------
             ! Get minimum boundary of the domain
             ! and set it into physics object.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) min_phys(1:num_dim)
             
             CALL physics_set_min_phys(phys,min_phys(1:num_dim),stat_info_sub)
             
             
          ELSE  IF (carg == 'MAX_PHYS') THEN
             
             !-----------------------------------------------
             ! Get maxmimum boundary of the domain.
             ! and set it into physics object.             
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) max_phys(1:num_dim)
             
             CALL physics_set_max_phys(phys,max_phys(1:num_dim),stat_info_sub)
             
             
          ELSE  IF (carg == 'LATTICE') THEN
             
             !-----------------------------------------------
             ! Get lattice type and
             ! and set it into physics object.             
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) lattice_type
             
             CALL physics_set_lattice_type(phys,lattice_type,stat_info_sub)
             
             
          ELSE  IF (carg == 'NUM_PART') THEN
             
             !-----------------------------------------------
             ! Get resolution, for particles 
             ! generated on lattice internally.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) num_part_dim(1:num_dim)
             
             CALL physics_set_num_part_dim(phys,&
                  num_part_dim(1:num_dim),stat_info_sub)
             
             
          ELSE  IF (carg == 'CUT_OFF') THEN
             
             !-----------------------------------------------
             ! Get cut off, i.e. compact support.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cut_off
             
             CALL physics_set_cut_off(phys,cut_off,stat_info_sub)
             
             
          ELSE  IF (carg == 'DT') THEN
           
             !-----------------------------------------------
             ! Get dt, if positive, 
             ! it will be used, otherwise,
             ! dt will be calculated internally
             ! according to certain creterion.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) dt
             
             CALL physics_set_dt(phys,dt,stat_info_sub)
             
             
          ELSE  IF (carg == 'STEP_START') THEN 
             
             !-----------------------------------------------
             ! Get step_start.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) step_start
             
             CALL physics_set_step_start(phys,step_start,stat_info_sub)
             
             
          ELSE  IF (carg == 'STEP_END') THEN 
             
             !-----------------------------------------------
             ! Get step_end.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) step_end
             
             CALL physics_set_step_end(phys,step_end,stat_info_sub)
             
             
          ELSE  IF (carg == 'TIME_START') THEN
             
             !-----------------------------------------------
             ! Get time_start.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) time_start
             
             CALL physics_set_time_start(phys,time_start,stat_info_sub)
             
             
          ELSE  IF (carg == 'TIME_END') THEN
             
             !-----------------------------------------------
             ! Get time_end.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) time_end
             
             CALL physics_set_time_end(phys,time_end,stat_info_sub)
             
             
             !#---------------------------------------------#
             !#---------------------------------------------#
             !    Get physics properties parameters.
             !#---------------------------------------------#             
             !#---------------------------------------------#
             
          ELSE IF (carg == 'RHO') THEN
             
             !-----------------------------------------------
             ! Mass density
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) rho
             
             CALL physics_set_rho(phys,rho,stat_info_sub)
             
             
          ELSE IF (carg =='ETA') THEN

             !-----------------------------------------------
             ! Shear viscosity
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) eta
             
             CALL physics_set_eta(phys,eta,stat_info_sub)
             
          ELSE IF (carg =='ETA_COEF') THEN

             !-----------------------------------------------
             ! Shear viscosity coefficient
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) eta_coef
             
             CALL physics_set_eta_coef(phys,eta_coef,stat_info_sub)
        
          ELSE IF (carg =='KSAI') THEN	  
             
             !-----------------------------------------------
             ! Bulk viscosity
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) ksai
             
             CALL physics_set_ksai(phys,ksai,stat_info_sub)
             
             
          ELSE IF (carg == 'KT') THEN
             
             !-----------------------------------------------
             ! KB (Boltzmann Constant) * 
             ! Temperature for SDPD.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) kt
             
             CALL physics_set_kt(phys,kt,stat_info_sub)
             
             
          ELSE IF (carg == 'C') THEN
             
             !-----------------------------------------------
             ! Sound speed
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) c
             
             CALL physics_set_c(phys,c,stat_info_sub)
             
             
          ELSE IF (carg == 'RHO_REF') THEN
             
             !-----------------------------------------------
             ! Reference mass density
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) rho_ref
             
             CALL physics_set_rho_ref(phys,rho_ref,stat_info_sub)
             
             
          ELSE IF (carg == 'GAMMA') THEN
             
             !-----------------------------------------------
             ! Power of density in state equation
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) gamma
             
             CALL physics_set_gamma(phys,gamma,stat_info_sub)
             
             
          ELSE IF (carg == 'BODY_FORCE_TYPE') THEN
             
             !----------------------------------------------
             ! Body force type
             !----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) body_force_type
             
             CALL physics_set_body_force_type(phys,&
                  body_force_type,stat_info_sub)
             
  
             !#---------------------------------------------#
             !#---------------------------------------------#
             !   Get relax run parameters.
             !#---------------------------------------------#             
             !#---------------------------------------------#
             
          ELSE IF (carg == 'RELAX_TYPE') THEN
             
             !-----------------------------------------------
             ! Relax run type
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) relax_type
             
             CALL physics_set_relax_type(phys,relax_type,stat_info_sub)

          ELSE  IF (carg == 'DT_RELAX') THEN
           
             READ(cvalue,*,IOSTAT=ios,ERR=200) dt_relax
             
             CALL physics_set_dt_relax(phys,dt,stat_info_sub)
             
          ELSE IF (carg == 'STEP_RELAX') THEN
             
             !-----------------------------------------------  
             ! Step for relax run.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) step_relax
             
             CALL physics_set_step_relax(phys,step_relax,stat_info_sub)
             

          ELSE IF (carg == 'TIME_RELAX') THEN
             
             !-----------------------------------------------
             ! Time for relax run.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) time_relax
             
             CALL physics_set_time_relax(phys,time_relax,stat_info_sub)

             
          ELSE IF (carg == 'DISORDER_LEVEL') THEN
             
             !-----------------------------------------------
             ! Time for relax run.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) disorder_level
             
             CALL physics_set_disorder_level(phys,disorder_level,stat_info_sub)


          ELSE IF (carg == 'KT_RELAX') THEN
             
             !-----------------------------------------------
             ! KB (Boltzmann Constant) * 
             ! Temperature for SDPD relax
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) kt_relax
             
             CALL physics_set_kt_relax(phys,kt_relax,stat_info_sub)

          ELSE IF (carg == 'C_RELAX') THEN
             
             !-----------------------------------------------
             ! Sound speed for SDPD relax
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) c_relax
             
             CALL physics_set_c_relax(phys,c_relax,stat_info_sub)
             
             
             !#---------------------------------------------#
             !#---------------------------------------------#
             !   Get viscoelastic parameters.
             !#---------------------------------------------#             
             !#---------------------------------------------#
           
             
          ELSE IF (carg == 'TAU_SM') THEN
             
             !-----------------------------------------------
             ! Relaxation time in stochastic model.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) tau_sm
             
             CALL physics_set_tau_sm(phys,tau_sm,stat_info_sub)

          ELSE IF (carg == 'K_SM') THEN
             
             !-----------------------------------------------
             ! Turbulence kinetic energy.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) k_sm
             
             CALL physics_set_k_sm(phys,k_sm,stat_info_sub)
             
           ELSE IF (carg == 'TAU') THEN
             
             !-----------------------------------------------
             ! Relaxation time of polymer molecules.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) tau
             
             CALL physics_set_tau(phys,tau,stat_info_sub)
            
          ELSE IF (carg == 'N_P') THEN
             
             !-----------------------------------------------
             ! number of molecules dumbbells.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) n_p
             
             CALL physics_set_n_p(phys, n_p,stat_info_sub)
             
             
          ELSE IF (carg == 'KT_P') THEN
             
             !-----------------------------------------------
             ! KB(Boltzmann Constant) * 
             ! Temperature for dumbell
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) kt_p
             
             CALL physics_set_kt_p(phys,kt_p,stat_info_sub)
             
             
          ELSE IF (carg == 'EIGEN_DYNAMICS') THEN
             
             !-----------------------------------------------
             ! Indicate if using eigen-dynamics.
             !-----------------------------------------------
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) eigen_dynamics
             
             CALL physics_set_eigen_dynamics(phys,&
                  eigen_dynamics,stat_info_sub)
             
             
          ELSE IF (carg == 'EVAL') THEN
             
             !-----------------------------------------------
             ! egenvalues for conformation tensor.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  eval(1:num_dim)
             
             CALL physics_set_eval(phys,eval(1:num_dim),stat_info_sub)
             
          ELSE IF (carg == 'EVAL_TOLERANCE') THEN
             
             !-----------------------------------------------
             ! egenvalues tolerance.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  eval_tolerance
             
             CALL physics_set_eval_tolerance(phys,eval_tolerance,stat_info_sub)

             
          ELSE IF (carg == 'EVEC') THEN
             
             !-----------------------------------------------
             ! egenvectors for conformation tensor.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                  evec(1:num_dim,1:num_dim)
             
             CALL physics_set_evec(phys,&
                  evec(1:num_dim,1:num_dim),stat_info_sub)

          ELSE IF (carg == 'EVEC_NORMALIZE') THEN
             
             !-----------------------------------------------
             ! Indicate if normalize eigenvectors.
             !-----------------------------------------------
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) evec_normalize
             
             CALL physics_set_evec_normalize(phys,&
                  evec_normalize,stat_info_sub)
             
          ELSE IF (carg == 'EVEC_TOLERANCE') THEN
             
             !-----------------------------------------------
             ! egenvector tolerance.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  evec_tolerance
             
             CALL physics_set_evec_tolerance(phys,evec_tolerance,stat_info_sub)

             !#---------------------------------------------#
             !#---------------------------------------------#
             !    Get body force and flow velocity.
             !#---------------------------------------------#             
             !#---------------------------------------------#
          
          ELSE IF (carg == 'BODY_FORCE') THEN
             
             !-----------------------------------------------
             ! Body force.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) body_force(1:num_dim)
             
             CALL physics_set_body_force(phys,&
                  body_force(1:num_dim),stat_info_sub)
             
             
          ELSE IF (carg == 'BODY_FORCE_D') THEN
             
             !-----------------------------------------------
             ! Increment of body force,
             ! in x-direction currently by default.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) body_force_d(1:num_dim)
             
             CALL physics_set_body_force_d(phys,&
                  body_force_d(1:num_dim),stat_info_sub)
             
             
          ELSE IF (carg == 'FLOW_DIRECTION') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) flow_direction
             
             CALL physics_set_flow_direction(phys,&
                  flow_direction,stat_info_sub)
             
             
          ELSE IF (carg == 'FLOW_WIDTH') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) flow_width
             
             CALL physics_set_flow_width(phys,flow_width,stat_info_sub)
             
             
          ELSE IF (carg == 'FLOW_V') THEN
             
             !-----------------------------------------------
             ! Flow velocity paramters, 
             ! in case we need fixed in-flow velocity.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) flow_v
             
             CALL physics_set_flow_v(phys,flow_v,stat_info_sub)
             
             
          ELSE IF (carg == 'FLOW_ADJUST_FREQ') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) flow_adjust_freq
             
             CALL physics_set_flow_adjust_freq(phys,&
                  flow_adjust_freq,stat_info_sub)
             
             !-----------------------------------------------
             ! Colloids  parameters 
             !-----------------------------------------------
             
          ELSE IF(carg == 'NUM_COLLOID' ) THEN
             
             !-----------------------------------------------
             ! Number of colloids.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) num_colloid

             CALL physics_set_num_colloid(phys,&
                  num_colloid,stat_info_sub)
          
             IF ( num_species > 1 .AND. num_colloid < 1 ) THEN
                PRINT *, "io_read_physics : ", &
                     "num species > 1, but num_colloid < 1"
                stat_info = -1
                GOTO 9999
             END IF
             
             coll_rho       = 0.0_MK
             coll_rho_type  = 0
             coll_translate = .TRUE.
             coll_rotate    = .TRUE.
             coll_place     = 1
             coll_noslip    = 2
             coll_body_force_type = 0.0
             coll_body_force(:) = 0.0_MK
             cc_lub_cut_off   = 0.0_MK
             cc_lub_cut_on    = 0.0_MK
             cc_repul_cut_off = 0.0_MK
             cc_repul_cut_on  = 0.0_MK
             cc_repul_F0      = 0.0_MK
             cw_lub_cut_off   = 0.0_MK
             cw_lub_cut_on    = 0.0_MK
             cw_repul_cut_off = 0.0_MK
             cw_repul_cut_on  = 0.0_MK
             cw_repul_F0      = 0.0_MK
             
             !-------------------------------------
             ! Allocate memory, assign zeros.
             !-------------------------------------
             
             ALLOCATE(coll_shape(num_colloid))
             ALLOCATE(coll_radius(num_dim,num_colloid))
             ALLOCATE(coll_freq(num_colloid))
             ALLOCATE(coll_m(num_colloid))
             ALLOCATE(coll_mmi(3,num_colloid))
             ALLOCATE(coll_x(num_dim,num_colloid))
             ALLOCATE(coll_v(num_dim,num_colloid))
             ALLOCATE(coll_acc_vector(4,num_colloid))
             ALLOCATE(coll_theta(3,num_colloid))
             ALLOCATE(coll_omega(3,num_colloid))
             
             coll_shape(:)     = 1
             coll_radius(:,:)  = 0.0_MK
             coll_freq(:)      = 0
             coll_m(:)         = 0.0_MK
             coll_mmi(:,:)     = 0.0_MK
             coll_x(:,:)       = 0.0_MK
             coll_v(:,:)       = 0.0_MK
             coll_acc_vector(:,:) = 0.0_MK
             coll_acc_vector(1,:) = 1.0_MK
             
             coll_theta(:,:)  = 0.0_MK
             coll_omega(:,:)  = 0.0_MK
             
             coll_index   = 0
             shape_index  = 0
             radius_index = 0
             freq_index   = 0
             m_index      = 0
             mmi_index    = 0
             x_index      = 0
             v_index      = 0
             rot_v_index    = 0
             theta_index    = 0
             omega_index  = 0
             
          ELSE IF (carg == 'COLL_RHO') THEN
             
             !-----------------------------------------------
             ! colloids mass density
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) coll_rho
             
          ELSE IF(carg == 'COLL_RHO_TYPE' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 ) THEN
             
             !-----------------------------------------------
             ! colloids density type
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) coll_rho_type
             
          ELSE IF (carg =='COLL_TRANSLATE' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 ) THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) coll_translate
             
          ELSE IF (carg =='COLL_ROTATE' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) coll_rotate
             
          ELSE IF(carg == 'COLL_PLACE' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) coll_place
        
          ELSE IF(carg == 'COLL_NOSLIP' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) coll_noslip
             
          ELSE IF (carg == 'COLL_BODY_FORCE_TYPE') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) coll_body_force_type
             
          ELSE IF (carg == 'COLL_BODY_FORCE' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 ) THEN
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                  coll_body_force(1:num_dim)
             
             !-----------------------------------------------
             !  colloid-colloid lubbrication cut off
             !-----------------------------------------------
             
          ELSE IF(carg == 'CC_LUB_CUT_OFF' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cc_lub_cut_off
             
             !-----------------------------------------------
             !  colloid-colloid lubrication cut on
             !-----------------------------------------------
             
          ELSE IF(carg == 'CC_LUB_CUT_ON' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cc_lub_cut_on
             
             !-----------------------------------------------
             !  colloid-colloid repulsive force cut off.
             !-----------------------------------------------
             
          ELSE IF(carg == 'CC_REPUL_CUT_OFF' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cc_repul_cut_off
             
             !-----------------------------------------------
             !  colloid-colloid repulsive force cut on
             !-----------------------------------------------
             
          ELSE IF(carg == 'CC_REPUL_CUT_ON' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cc_repul_cut_on
             
             !-----------------------------------------------
             !  colloid-colloid maximum repulsive force.
             !-----------------------------------------------
             
          ELSE IF(carg == 'CC_REPUL_F0' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cc_repul_F0
             
             !-----------------------------------------------
             ! colloid-wall lubrication cut off
             !-----------------------------------------------
             
          ELSE IF(carg == 'CW_LUB_CUT_OFF' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cw_lub_cut_off
             
             !-----------------------------------------------
             ! colloid-wall lubrication cut on
             !-----------------------------------------------
             
          ELSE IF(carg == 'CW_LUB_CUT_ON' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cw_lub_cut_on
             
             !-----------------------------------------------
             ! colloid-wall repulsive force cut off.
             !-----------------------------------------------
             
          ELSE IF(carg == 'CW_REPUL_CUT_OFF' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cw_repul_cut_off
             
             !-----------------------------------------------
             ! colloid-wall repulsive force cut on.
             !-----------------------------------------------
             
          ELSE IF(carg == 'CW_REPUL_CUT_ON' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cw_repul_cut_on
             
             !-----------------------------------------------
             ! colloid-wall repulsive force maximum
             !-----------------------------------------------
             
          ELSE IF(carg == 'CW_REPUL_F0' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0) THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cw_repul_F0
             
             !-----------------------------------------------
             ! Individual properties of colloids start
             !-----------------------------------------------
             
          ELSE IF (carg == 'COLL_SHAPE' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 .AND. &
               shape_index >= 0 .AND.  &
               shape_index < num_colloid ) THEN
             
             shape_index = shape_index + 1
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) &
                  coll_shape(shape_index)
             
          ELSE IF (carg == 'COLL_RADIUS' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 .AND. &
               radius_index >= 0 .AND.  &
               radius_index <num_colloid ) THEN
             
             radius_index = radius_index + 1
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                  coll_radius(1:num_dim,radius_index)
             
          ELSE IF (carg == 'COLL_FREQ' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 .AND. &
               freq_index >= 0 .AND.  &
               freq_index <num_colloid ) THEN
             
             freq_index = freq_index + 1
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                  coll_freq(freq_index)
             
          ELSE IF (carg == 'COLL_M' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 .AND. &
               m_index >= 0 .AND.  &
               m_index <num_colloid ) THEN
             
             m_index = m_index + 1
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                  coll_m(m_index)
             
          ELSE IF (carg == 'COLL_MMI' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 .AND. &
               mmi_index >= 0 .AND.  &
               mmi_index <num_colloid ) THEN
             
             mmi_index = mmi_index + 1
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                  coll_mmi(1:3,mmi_index)
             
          ELSE IF (carg == 'COLL_X' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 .AND. &
               x_index >= 0 .AND.  &
               x_index <num_colloid ) THEN
             
             x_index = x_index + 1

             READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                  coll_x(1:num_dim,x_index)
             
          ELSE IF (carg == 'COLL_V' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 .AND. &
               v_index >= 0 .AND.  &
               v_index <num_colloid ) THEN
             
             v_index = v_index + 1
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) &
                  coll_v(1:num_dim,v_index)

          ELSE IF (carg == 'COLL_ROT_VECTOR' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 .AND. &
               rot_v_index >= 0 .AND.  &
               rot_v_index <num_colloid ) THEN
             
             rot_v_index = rot_v_index + 1
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                  coll_acc_vector(1:4,rot_v_index)
       
          ELSE IF (carg == 'COLL_THETA' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 .AND. &
               theta_index >= 0 .AND.  &
               theta_index <num_colloid ) THEN
             
             theta_index = theta_index + 1
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                  coll_theta(1:3,theta_index)
             
          ELSE IF (carg == 'COLL_OMEGA' .AND. &
               num_species > 1 .AND. &
               num_colloid > 0 .AND. &
               omega_index >= 0 .AND.  &
               omega_index <num_colloid ) THEN
             
             omega_index = omega_index + 1
             
             READ(cvalue,*,IOSTAT=ios, ERR=200) &
                  coll_omega(1:3,omega_index)
             
             !#---------------------------------------------#
             !#---------------------------------------------#
             ! Get boundary parameters.
             !#---------------------------------------------#             
             !#---------------------------------------------#
             
          ELSE IF (carg == 'BCDEF') THEN
             
             !-----------------------------------------------
             ! Physics boundary conditions.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  bcdef(1:2*num_dim)
             
             CALL physics_set_bcdef(phys,bcdef(1:2*num_dim),stat_info_sub)
             
             
          ELSE IF (carg == 'SHEAR_TYPE') THEN
             
             !-----------------------------------------------
             ! shear boundary type.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  shear_type(1:2*num_dim)
          

          ELSE IF (carg == 'SHEAR_V') THEN
             
             !-----------------------------------------------
             ! shear boundary's initial rate
             !
             ! 2D : shear velocity of x side, in y direction;
             !      shear velocity of y side, in x direction;
             !
             !-----------------------------------------------
             
             IF ( num_dim == 2 ) THEN
                
                READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                     shear_v(2,1),shear_v(2,2), &
                     shear_v(1,3),shear_v(1,4)
                
             ELSE IF ( num_dim == 3 ) THEN
                
                READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                     shear_v(2,1),shear_v(3,1), &
                     shear_v(2,2),shear_v(3,2), &
                     shear_v(1,3),shear_v(3,3), &
                     shear_v(1,4),shear_v(3,4), &
                     shear_v(1,5),shear_v(2,5), &
                     shear_v(1,6),shear_v(2,6)
                 
             END IF
             
             
          ELSE IF (carg == 'SHEAR_FREQ') THEN
             
             !-----------------------------------------------
             ! shear boundary's
             ! oscillating frequncy
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  shear_freq(1:2*num_dim)
             

          ELSE IF (carg == 'WALL_RHO_TYPE') THEN
             
             !-----------------------------------------------
             ! Wall boundary density type.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  wall_rho_type
             

          ELSE IF (carg == 'WALL_NOSLIP') THEN
             
             !-----------------------------------------------
             ! Wall boundary no slip type.
             !-----------------------------------------------
             
             READ(cvalue,*,IOSTAT=ios, ERR=200)  wall_noslip
             
          END IF ! check argument name
          
        END DO ! read each line
       
        !----------------------------------------------------
        ! End of file.
        !----------------------------------------------------
        
        
200     CONTINUE
        
        !----------------------------------------------------
        ! Something went wrong.
        !----------------------------------------------------
        
        WRITE(cbuf,'(A,I5,2A)') 'Error reading line ',iline,     &
             &  ' of file : ',this%physics_config_file(1:ilenphys)
        
        ilen = LEN_TRIM(cbuf)
        PRINT *,'io_read_physics_config_file : ', cbuf(1:ilen)
        
        stat_info = -1
        GOTO 9999
        
        
9999	CONTINUE
	
        
        IF ( stat_info == 0 ) THEN
           
           !-------------------------------------------------
           ! If everything went smooth,
           ! create boundary object.
           !-------------------------------------------------
           
           CALL boundary_new(tboundary,&
                num_dim,stat_info_sub)
           CALL boundary_set_bcdef(tboundary,&
                bcdef(1:2*num_dim),stat_info_sub)
           CALL boundary_set_shear_type(tboundary,&
                shear_type(1:2*num_dim),stat_info_sub)
           CALL boundary_set_shear_v0(tboundary,&
                shear_v(1:num_dim,1:2*num_dim),stat_info_sub)
           CALL boundary_set_shear_freq(tboundary,&
                shear_freq(1:2*num_dim),stat_info_sub)
           CALL boundary_set_wall_rho_type(tboundary, &
                wall_rho_type, stat_info_sub)
           CALL boundary_set_wall_noslip_type(tboundary,&
                wall_noslip,stat_info_sub)
           CALL boundary_set_min_phys(tboundary,&
                min_phys(1:num_dim),stat_info_sub)
           CALL boundary_set_max_phys(tboundary,&
                max_phys(1:num_dim),stat_info_sub)
           
           !-------------------------------------------------
           ! shear_rate = (shear_v(min) - shear_v(max))/L
           !-------------------------------------------------
           
           DO i = 1, num_dim
              
              shear_rate(1:num_dim,i) = &
                   ( shear_v(1:num_dim,2*i-1) -  &
                   shear_v(1:num_dim,2*i) ) / &
                   ( max_phys(i) - min_phys(i) )
              
           END DO
           
           CALL boundary_set_shear_rate(tboundary,&
                shear_rate(1:num_dim,1:num_dim),stat_info_sub)
           CALL physics_set_boundary(phys,&
                tboundary,stat_info_sub)
           
           IF ( num_colloid > 0 ) THEN
              
              !----------------------------------------------
              ! Set properties of colloids.
              !----------------------------------------------
              
              CALL colloid_new(colloids,&
                   num_dim,num_colloid,stat_info_sub)
              CALL colloid_set_rho(colloids, &
                   coll_rho,stat_info_sub)         
              CALL colloid_set_rho_type(colloids, &
                   coll_rho_type,stat_info_sub)
              CALL colloid_set_translate(colloids,&
                   coll_translate,stat_info_sub)
              CALL colloid_set_rotate(colloids,&
                   coll_rotate,stat_info_sub) 
              CALL colloid_set_place(colloids,&
                   coll_place,stat_info_sub)
              CALL colloid_set_noslip_type(colloids,&
                   coll_noslip,stat_info_sub)
              CALL colloid_set_body_force_type(colloids,&
                   coll_body_force_type,stat_info_sub)             
              CALL colloid_set_body_force(colloids,&
                   coll_body_force(1:num_dim),&
                   stat_info_sub)
             
              CALL colloid_set_cc_lub_cut_off(colloids, &
                   cc_lub_cut_off,stat_info_sub)
              CALL colloid_set_cc_lub_cut_on(colloids, &
                   cc_lub_cut_on,stat_info_sub)
              CALL colloid_set_cc_repul_cut_off(colloids, &
                   cc_repul_cut_off,stat_info_sub)
              CALL colloid_set_cc_repul_cut_on(colloids, &
                   cc_repul_cut_on,stat_info_sub)            
              CALL colloid_set_cc_repul_F0(colloids, &
                   cc_repul_F0,stat_info_sub)
              
              CALL colloid_set_cw_lub_cut_off(colloids, &
                   cw_lub_cut_off,stat_info_sub)
              CALL colloid_set_cw_lub_cut_on(colloids, &
                   cw_lub_cut_on,stat_info_sub)
              CALL colloid_set_cw_repul_cut_off(colloids, &
                   cw_repul_cut_off,stat_info_sub)
              CALL colloid_set_cw_repul_cut_on(colloids, &
                   cw_repul_cut_on,stat_info_sub)            
              CALL colloid_set_cw_repul_F0(colloids, &
                   cw_repul_F0,stat_info_sub)
              
              CALL colloid_set_shape(colloids,&
                   coll_shape(1:num_colloid),stat_info_sub)
              CALL colloid_set_radius(colloids,&
                   coll_radius(1:num_dim,1:num_colloid),stat_info_sub)
              CALL colloid_set_freq(colloids,&
                   coll_freq(1:num_colloid),stat_info_sub)
              CALL colloid_set_m(colloids,&
                   coll_m(1:num_colloid),stat_info_sub)
              CALL colloid_set_mmi(colloids,&
                   coll_mmi(1:3,1:num_colloid),stat_info_sub)
              CALL colloid_set_x(colloids,&
                   coll_x(1:num_dim,1:num_colloid),stat_info_sub)
              CALL colloid_set_v(colloids,&
                   coll_v(1:num_dim,1:num_colloid),stat_info_sub)
              
              !CALL colloid_set_rotation_vector(colloids,&
              !     coll_acc_vector(1:4,1:num_colloid),stat_info_sub)           
              CALL colloid_set_accumulation_vector(colloids,&
                   coll_acc_vector(1:4,1:num_colloid),stat_info_sub)           
           
                    
              CALL colloid_set_theta(colloids,&
                   coll_theta(1:3,1:num_colloid),stat_info_sub)           
              CALL colloid_set_omega(colloids,&
                   coll_omega(1:3,1:num_colloid),stat_info_sub)
              
              !----------------------------------------------
              ! Record some often used physics properties
              ! in Colloid class object.
              !----------------------------------------------
              
              CALL colloid_set_min_phys(colloids,&
                   min_phys(1:num_dim),stat_info_sub)
              CALL colloid_set_max_phys(colloids,&
                   max_phys(1:num_dim),stat_info_sub)
              CALL colloid_set_bcdef(colloids,&
                   bcdef(1:2*num_dim),stat_info_sub)
              CALL colloid_set_boundary(colloids,tboundary,&
                   stat_info_sub)
              CALL colloid_set_cut_off(colloids,cut_off, &
                   stat_info_sub)
              CALL colloid_set_eta(colloids, &
                   eta,stat_info_sub)
              
              CALL physics_set_colloid(phys,&
                   colloids,stat_info_sub)
              
              
           END IF
           
        END IF
        
	!----------------------------------------------------
      	! Release dynamic memory and Close file.
      	!----------------------------------------------------
        
        IF(ALLOCATED(coll_shape))THEN
           DEALLOCATE(coll_shape)
        END IF
        
        IF(ALLOCATED(coll_radius))THEN
           DEALLOCATE(coll_radius)
        END IF

        IF(ALLOCATED(coll_freq))THEN
           DEALLOCATE(coll_freq)
        END IF
        
        IF(ALLOCATED(coll_m))THEN
           DEALLOCATE(coll_m)
        END IF

        IF(ALLOCATED(coll_mmi))THEN
           DEALLOCATE(coll_mmi)
        END IF
                
        IF(ALLOCATED(coll_x))THEN
           DEALLOCATE(coll_x)
        END IF
        
        IF(ALLOCATED(coll_v))THEN
           DEALLOCATE(coll_v)
        END IF
        
        IF(ALLOCATED(coll_acc_vector))THEN
           DEALLOCATE(coll_acc_vector)
        END IF

        IF(ALLOCATED(coll_theta))THEN
           DEALLOCATE(coll_theta)
        END IF
        
        IF(ALLOCATED(coll_omega))THEN
           DEALLOCATE(coll_omega)
        END IF
        
        
        CLOSE(this%physics_config_unit)
        
        
#ifdef __DEBUG
        IF ( debug_threshold < 3) THEN
           PRINT *, 'io_read_physics_config : ended'
        END IF
#endif
        
        RETURN        
        
      END SUBROUTINE io_read_physics_config
      
      
      
