      SUBROUTINE io_write_restart_physics(this,&
           rank,step,time,d_physics,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_write_restart_physics
        !----------------------------------------------------
        !
        ! Purpose     : Write essential physics parameters
        !               (including colloid) to be able to 
        !               restart the simulation again.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.2 18.11.2010, some of physics
        !               parameters are not needed for
        !               writing restart, e.g., 
        !               for non-Newtonian fluid, all the
        !               viscoelastic parameters are not
        !               needed.
        !
        !               V0.1 01.04.2009, original version.
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
	
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)             :: rank
        INTEGER,  INTENT(IN)	        :: step
        REAL(MK), INTENT(IN)	        :: time
        TYPE(Physics), INTENT(IN)       :: d_physics
        INTEGER,  INTENT(OUT)	        :: stat_info
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub

        !----------------------------------------------------
        ! Control parameters.
        !----------------------------------------------------
        
        LOGICAL                         :: Newtonian
        
        !----------------------------------------------------
        ! physics parameters.
        !----------------------------------------------------
        
        INTEGER                         :: num_species
        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(:), POINTER :: min_phys
        REAL(MK), DIMENSION(:), POINTER :: max_phys
        INTEGER, DIMENSION(:), POINTER  :: num_part_dim
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
        REAL(MK)                        :: n_p
        REAL(MK)                        :: kt_p
        LOGICAL                         :: eigen_dynamics
        REAL(MK),DIMENSION(:),POINTER   :: eval
        REAL(MK)                        :: eval_tolerance
        REAL(MK),DIMENSION(:,:),POINTER :: evec
        LOGICAL                         :: evec_normalize
        REAL(MK)                        :: evec_tolerance

        INTEGER                         :: body_force_type
        REAL(MK), DIMENSION(:), POINTER :: body_force
        REAL(MK), DIMENSION(:), POINTER :: body_force_d
        INTEGER                         :: flow_direction
        REAL(MK)                        :: flow_width
        REAL(MK)                        :: flow_v
        INTEGER                         :: flow_adjust_freq

        !----------------------------------------------------
        ! Colloid parameters.
        !----------------------------------------------------

        INTEGER                         :: num_colloid
        TYPE(Colloid),POINTER           :: colloids
        INTEGER                         :: coll_rho_type
        LOGICAL                         :: coll_translate
        LOGICAL                         :: coll_rotate
        INTEGER                         :: coll_place
        INTEGER                         :: coll_noslip
        INTEGER                         :: coll_body_force_type
        REAL(MK), DIMENSION(:), POINTER :: coll_body_force
        REAL(MK)                        :: coll_cc_lub_cut_off
        REAL(MK)                        :: coll_cc_lub_cut_on
        REAL(MK)                        :: coll_cc_repul_cut_off
        REAL(MK)                        :: coll_cc_repul_cut_on
        REAL(MK)                        :: coll_cc_repul_F0
        REAL(MK)                        :: coll_cw_lub_cut_off
        REAL(MK)                        :: coll_cw_lub_cut_on    
        REAL(MK)                        :: coll_cw_repul_cut_off
        REAL(MK)                        :: coll_cw_repul_cut_on
        REAL(MK)                        :: coll_cw_repul_F0
       
        INTEGER, DIMENSION(:), POINTER  :: coll_shape
        REAL(MK),DIMENSION(:,:),POINTER :: coll_radius
        INTEGER, DIMENSION(:), POINTER  :: coll_freq
        REAL(MK), DIMENSION(:), POINTER :: coll_m
        REAL(MK),DIMENSION(:,:),POINTER :: coll_mmi
        REAL(MK),DIMENSION(:,:),POINTER :: coll_x
        REAL(MK),DIMENSION(:,:),POINTER :: coll_v
        REAL(MK),DIMENSION(:,:),POINTER :: coll_acc_vector
        REAL(MK),DIMENSION(:,:),POINTER :: coll_theta        
        REAL(MK),DIMENSION(:,:),POINTER :: coll_omega
        
        !----------------------------------------------------
        ! Boundary parameters.
        !----------------------------------------------------
        
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER, DIMENSION(:), POINTER  :: shear_type
        REAL(MK),DIMENSION(:,:),POINTER :: shear_v
        REAL(MK),DIMENSION(:),POINTER   :: shear_freq
        INTEGER                         :: wall_rho_type
        INTEGER                         :: wall_noslip

        LOGICAL                         :: lexist
        CHARACTER(len=MAX_CHAR)	        :: file_name
        CHARACTER(len=2*MAX_CHAR)	:: cbuf
        INTEGER			        :: i,j

        !----------------------------------------------------
        ! Initialization.
        !----------------------------------------------------
        
        IF ( rank /= 0 ) THEN
           PRINT *, "io_write_restart_physics : ", &
                "Can only be called by root processor !"
           stat_info = -1
           GOTO 9999
        END IF
        
        stat_info     = 0
        stat_info_sub = 0
        
        Newtonian = &
             control_get_Newtonian(this%ctrl,stat_info_sub)
        
        NULLIFY(min_phys)
        NULLIFY(max_phys)
        NULLIFY(num_part_dim)
        NULLIFY(body_force)
        NULLIFY(body_force_d)
        
        NULLIFY(colloids)
        NULLIFY(coll_body_force)
        coll_translate = .FALSE.
        coll_rotate    = .FALSE.
        coll_place     = 1
        coll_noslip    = 2
        coll_cc_lub_cut_off = 0.0_MK
        coll_cc_lub_cut_on  = 0.0_MK
        coll_cc_repul_cut_off = 0.0_MK
        coll_cc_repul_cut_on  = 0.0_MK
        coll_cc_repul_F0      = 0.0_MK
        coll_cw_lub_cut_off = 0.0_MK
        coll_cw_lub_cut_on  = 0.0_MK
        coll_cw_repul_cut_off = 0.0_MK
        coll_cw_repul_cut_on  = 0.0_MK
        coll_cw_repul_F0      = 0.0_MK

        NULLIFY(coll_shape)
        NULLIFY(coll_radius)
        NULLIFY(coll_freq)
        NULLIFY(coll_m)
        NULLIFY(coll_mmi)
        NULLIFY(coll_x)
        NULLIFY(coll_v)
        NULLIFY(coll_acc_vector)
        NULLIFY(coll_theta)
        NULLIFY(coll_omega)   
        
        
        NULLIFY(bcdef)
        NULLIFY(tboundary)
        NULLIFY(shear_v)
        NULLIFY(shear_type)
        NULLIFY(shear_freq)
        CALL physics_get_boundary(d_physics,tboundary,stat_info_sub)
        CALL boundary_get_shear_type(tboundary,shear_type,stat_info_sub)
        CALL boundary_get_shear_v(tboundary,shear_v,stat_info_sub)
        CALL boundary_get_shear_freq(tboundary,shear_freq,stat_info_sub)
        wall_rho_type = &
             boundary_get_wall_rho_type(tboundary,stat_info_sub)
        walL_noslip   = &
             boundary_get_wall_noslip_type(tboundary,stat_info_sub)
        
        WRITE(file_name,'(A,I8.8,A)') &
             TRIM(this%restart_physics_file),step,'.mcf'
        
        INQUIRE(FILE=file_name,EXIST=lexist)
        
        IF( lexist ) THEN
           
           OPEN(UNIT=this%restart_physics_unit,FILE=file_name,&
                STATUS="REPLACE",FORM=this%restart_physics_fmt,&
                ACTION="WRITE",  POSITION="APPEND",IOSTAT=stat_info_sub)
           PRINT *, "io_write_restart_physics : ", &
                "Existing restart_physics file is replaced !"
           
        ELSE
           
           OPEN(UNIT=this%restart_physics_unit,FILE=file_name,&
                STATUS="NEW",FORM=this%restart_physics_fmt,&
                ACTION="WRITE", IOSTAT=stat_info_sub)
           
        END IF
        
        
        IF(stat_info_sub /= 0) THEN
           PRINT *, "io_write_restart_physics : ",&
                "creating/opening the ", &
                TRIM(this%restart_physics_file), " has error!"
           stat_info = -1
           GOTO 9999           
        END IF
        
        !------------------------------------------------------------
        ! Write Comments in the restart physics file.
        !------------------------------------------------------------
        
        WRITE(cbuf, '(A)') '#---------------------------------------------'
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(cbuf, '(A)') '# The meaning of following variables can be '
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(cbuf, '(A)') '# looked up in original physics_config.mcf file'
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(cbuf, '(A)') '#---------------------------------------------'
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '

        !-------------------------------------------------
        ! Wirte number of species
        !-------------------------------------------------
        
        num_species = physics_get_num_species(d_physics,stat_info_sub)
        WRITE(cbuf, '(A,I2)') 'num_species = ', num_species
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '

        !-------------------------------------------------
        ! Wirte number of dimensionality
        !-------------------------------------------------
        
        num_dim = physics_get_num_dim(d_physics,stat_info_sub)
        WRITE(cbuf, '(A,I2)') 'num_dim = ', num_dim         
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        !-------------------------------------------------
        ! Wirte minimum physiscs of the domain
        !-------------------------------------------------
        
        CALL physics_get_min_phys(d_physics,min_phys,stat_info_sub)
        IF(num_dim == 2) THEN
           WRITE(cbuf, '(2(A,E16.8))') &
                'min_phys = ',min_phys(1),',', min_phys(2)
        ELSE IF (num_dim == 3 ) THEN
           WRITE(cbuf, '(3(A,E16.8))') &
                'min_phys = ',min_phys(1),&
                ',', min_phys(2), ',', min_phys(3)
        END IF
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        !-------------------------------------------------
        ! Wirte maximum physiscs of the domain.
        !-------------------------------------------------
        
        CALL physics_get_max_phys(d_physics,max_phys,stat_info_sub)
        IF(num_dim == 2) THEN
           WRITE(cbuf, '(2(A,E16.8))') &
                'max_phys = ',max_phys(1),',', max_phys(2)
        ELSE IF (num_dim == 3 ) THEN
           WRITE(cbuf, '(3(A,E16.8))') &
                'max_phys = ',max_phys(1),&
                ',', max_phys(2), ',', max_phys(3)
        END IF
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        !-------------------------------------------------
        ! Write resolution 
        ! (this doesn't make much sence, since partciles' 
        ! positions are chaotic after simulation started,
        ! not on lattice any more)
        !-------------------------------------------------
        
        CALL physics_get_num_part_dim(d_physics,num_part_dim,stat_info_sub)
        
        IF(num_dim == 2) THEN
           WRITE(cbuf, '(2(A,I4))') &
                'num_part = ',num_part_dim(1),',', num_part_dim(2)
        ELSE IF (num_dim == 3 ) THEN
           WRITE(cbuf, '(3(A,I4))') &
                'num_part = ',num_part_dim(1),&
                ',', num_part_dim(2), ',', num_part_dim(3)
        END IF
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        
        !----------------------------------------------------
        ! Write cut off.
        !-----------------------------------------------------
        
        cut_off = physics_get_cut_off(d_physics,stat_info_sub)
        WRITE(cbuf, '(A,E16.8)') 'cut_off = ', cut_off
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
   
        
        !----------------------------------------------------
        ! Write  dt.
        !----------------------------------------------------
        
        dt = physics_get_dt(d_physics,stat_info_sub)
        WRITE(cbuf, '(A,E16.8)') 'dt = ', dt
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        !----------------------------------------------------
        ! Write step_start
        !-----------------------------------------------------
        
        step_start = step
        WRITE(cbuf, '(A,I10)') 'step_start = ', step_start
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
   
        !----------------------------------------------------
        ! Write step_end
        !-----------------------------------------------------
          
        step_end = physics_get_step_end(d_physics,stat_info_sub)
        WRITE(cbuf, '(A,I10)') 'step_end = ', step_end
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
   
        !----------------------------------------------------
        ! Write time_start
        !-----------------------------------------------------
        
        time_start = time
        WRITE(cbuf, '(A,E16.8)') 'time_start = ', time_start 
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        !----------------------------------------------------
        ! Write time_end
        !-----------------------------------------------------
        
        time_end = physics_get_time_end(d_physics,stat_info_sub)
        WRITE(cbuf, '(A,E16.8)') 'time_end = ', time_end
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        !----------------------------------------------------
        ! Write physics input parameters 
        !
        ! rho      : mass density
        ! eta      : dynamic viscosity
        ! ksai     : bulk viscosity
        ! kt       : Boltzmann constant * temperature for SDPD.
        ! c        : sound speed
        ! rho_ref  : reference density
        ! gamma    : power in state equation
        !----------------------------------------------------
        
        rho = physics_get_rho(d_physics,stat_info_sub)        
        WRITE(cbuf, '(A,E16.8)') 'rho = ', rho
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        eta = physics_get_eta(d_physics,stat_info_sub)        
        WRITE(cbuf, '(A,E16.8)') 'eta = ', eta
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        eta_coef = physics_get_eta_coef(d_physics,stat_info_sub)        
        WRITE(cbuf, '(A,E16.8)') 'eta_coef = ', eta_coef
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
     
        ksai = physics_get_ksai(d_physics,stat_info_sub)        
        WRITE(cbuf, '(A,E16.8)') 'ksai = ', ksai
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
   
        kt = physics_get_kt(d_physics,stat_info_sub) 
        WRITE(cbuf, '(A,E16.8)') 'kt = ', kt
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        c = physics_get_c(d_physics,stat_info_sub) 
        WRITE(cbuf, '(A,E16.8)') 'c = ', c
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        rho_ref = physics_get_rho_ref(d_physics,stat_info_sub)        
        WRITE(cbuf, '(A,E16.8)') 'rho_ref = ', rho_ref
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        gamma = physics_get_gamma(d_physics,stat_info_sub)        
        WRITE(cbuf, '(A,E16.8)') 'gamma = ', gamma
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '

        
        !--------------------------------------------------
        ! Write relax parameters.
        !--------------------------------------------------
        
        relax_type = physics_get_relax_type(d_physics,stat_info_sub)
        WRITE(cbuf, '(A,I10)') 'relax_type = ', relax_type
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        dt_relax = physics_get_dt_relax(d_physics,stat_info_sub)
        WRITE(cbuf, '(A,E16.8)') 'dt_relax = ', dt_relax
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
     
        step_relax = physics_get_step_relax(d_physics,stat_info_sub)
        WRITE(cbuf, '(A,I10)') 'step_relax = ', step_relax
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        time_relax = physics_get_time_relax(d_physics,stat_info_sub) 
        WRITE(cbuf, '(A,E16.8)') 'time_relax = ', time_relax
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        disorder_level = physics_get_disorder_level(d_physics,stat_info_sub) 
        WRITE(cbuf, '(A,E16.8)') 'disorder_level = ', disorder_level
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
     
        kt_relax = physics_get_kt_relax(d_physics,stat_info_sub) 
        WRITE(cbuf, '(A,E16.8)') 'kt_relax = ', kt_relax
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '

        c_relax = physics_get_c_relax(d_physics,stat_info_sub) 
        WRITE(cbuf, '(A,E16.8)') 'c_relax = ', c_relax
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
     
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '

        
        !----------------------------------------------------
        ! Write egenvector dynamics for Non-Newtonian fluids
        ! oldroyd-B model.
        !----------------------------------------------------
        
        IF ( .NOT. Newtonian ) THEN
           
           tau = physics_get_tau(d_physics,stat_info_sub)        
           WRITE(cbuf, '(A,E16.8)') 'tau = ', tau
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           n_p = physics_get_n_p(d_physics,stat_info_sub)        
           WRITE(cbuf, '(A,E16.8)') 'n_p = ', n_p
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           kt_p = physics_get_kt_p(d_physics,stat_info_sub) 
           WRITE(cbuf, '(A,E16.8)') 'kt_p = ', kt_p
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           eigen_dynamics = physics_get_eigen_dynamics(d_physics,stat_info_sub)
           WRITE(cbuf, '(A,L)') 'eigen_dynamics = ', eigen_dynamics
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
           CALL physics_get_eval(d_physics,eval,stat_info_sub)
           
           IF(num_dim == 2) THEN
              WRITE(cbuf, '(2(A,E16.8))') &
                   'eval = ',eval(1),',', eval(2)
           ELSE IF (num_dim == 3 ) THEN
              WRITE(cbuf, '(3(A,E16.8))') &
                   'eval = ',eval(1),',', eval(2),',', eval(3)
           END IF
           
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           eval_tolerance = physics_get_eval_tolerance(d_physics,stat_info_sub)        
           WRITE(cbuf, '(A,E16.8)') 'eval_tolerance = ', eval_tolerance
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
           CALL physics_get_evec(d_physics,evec,stat_info_sub)
           
           IF(num_dim == 2) THEN
              WRITE(cbuf, '(4(A,E16.8))') &
                   'evec = ',evec(1,1),',', evec(2,1),&
                   ',', evec(1,2),',', evec(2,2)
           ELSE IF (num_dim == 3 ) THEN
              WRITE(cbuf, '(9(A,E16.8))') &
                   'evec = ',evec(1,1),',', evec(2,1),',', evec(3,1),&
                   ( (',', evec(i,j),i=1,3),j=2,3 )
           END IF
           
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           evec_normalize = physics_get_evec_normalize(d_physics,stat_info_sub)
           WRITE(cbuf, '(A,L)') 'evec_normalize = ', evec_normalize
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
           evec_tolerance = physics_get_evec_tolerance(d_physics,stat_info_sub)        
           WRITE(cbuf, '(A,E16.8)') 'evec_tolerance = ', evec_tolerance
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
        END IF  ! non-Newtonian
        
        !--------------------------------------------------
        ! Write external body force
        !--------------------------------------------------
        
        body_force_type = physics_get_body_force_type(d_physics,stat_info_sub)
        
        WRITE(cbuf, '(A,I3)') 'body_force_type = ', body_force_type
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        CALL physics_get_body_force(d_physics,body_force,stat_info_sub)
        
        IF(num_dim == 2) THEN
           WRITE(cbuf, '(2(A,E16.8))') &
                'body_force = ',body_force(1),',', body_force(2)
        ELSE IF (num_dim == 3 ) THEN
           WRITE(cbuf, '(3(A,E16.8))') &
                'body_force = ',body_force(1),&
                ',', body_force(2), ',', body_force(3)
        END IF
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        CALL physics_get_body_force_d(d_physics,body_force_d,stat_info_sub)
        
        IF(num_dim == 2) THEN
           WRITE(cbuf, '(2(A,E16.8))') &
                'body_force_d = ',body_force_d(1),',', body_force_d(2)
        ELSE IF (num_dim == 3 ) THEN
           WRITE(cbuf, '(3(A,E16.8))') &
                'body_force_d = ',body_force_d(1),&
                ',', body_force_d(2), ',', body_force_d(3)
        END IF
        
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        !--------------------------------------------------
        ! Write flow velocity parameters
        !--------------------------------------------------
        
        flow_direction = physics_get_flow_direction(d_physics,stat_info_sub)        
        WRITE(cbuf, '(A,I5)') 'flow_direction = ', flow_direction
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        flow_width = physics_get_flow_width(d_physics,stat_info_sub)        
        WRITE(cbuf, '(A,E16.8)') 'flow_width = ', flow_width
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        flow_v = physics_get_flow_v(d_physics,stat_info_sub)        
        WRITE(cbuf, '(A,E16.8)') 'flow_v = ', flow_v
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
    
        flow_adjust_freq = physics_get_flow_adjust_freq(d_physics,stat_info_sub)        
        WRITE(cbuf, '(A,I5)') 'flow_adjust_freq = ', flow_adjust_freq
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        
        !----------------------------------------------------
        ! Wirte colloids  parameters.
        !----------------------------------------------------
        
        num_colloid = physics_get_num_colloid(d_physics,stat_info_sub) 
        
        WRITE(cbuf, '(A,I3)') 'num_colloid = ', num_colloid
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        IF( num_species > 1 .AND. num_colloid > 0 ) THEN 
        
           CALL physics_get_colloid(d_physics,colloids,stat_info_sub)
           
           coll_rho_type  = &
                colloid_get_rho_type(colloids,stat_info_sub)
           coll_translate = &
                colloid_get_translate(colloids,stat_info_sub)
           coll_rotate    = &
                colloid_get_rotate(colloids,stat_info_sub)
           coll_place     = &
                colloid_get_place(colloids,stat_info_sub)
           coll_noslip    = &
                colloid_get_noslip_type(colloids,stat_info_sub)
           coll_body_force_type= &
                colloid_get_body_force_type(colloids,stat_info_sub)
           CALL colloid_get_body_force(colloids,coll_body_force,stat_info_sub)
           coll_cc_lub_cut_off = &
                colloid_get_cc_lub_cut_off(colloids,stat_info_sub)
           coll_cc_lub_cut_on = &
                colloid_get_cc_lub_cut_on(colloids,stat_info_sub)
           coll_cc_repul_cut_off = &
                colloid_get_cc_repul_cut_off(colloids,stat_info_sub)
           coll_cc_repul_cut_on = &
                colloid_get_cc_repul_cut_on(colloids,stat_info_sub)         
           coll_cc_repul_F0      = &
                colloid_get_cc_repul_F0(colloids,stat_info_sub)  
           coll_cw_lub_cut_off = &
                colloid_get_cw_lub_cut_off(colloids,stat_info_sub)
           coll_cw_lub_cut_on = &
                colloid_get_cw_lub_cut_on(colloids,stat_info_sub)         
           coll_cw_repul_cut_off = &
                colloid_get_cw_repul_cut_off(colloids,stat_info_sub)
           coll_cw_repul_cut_on = &
                colloid_get_cw_repul_cut_on(colloids,stat_info_sub)        
           coll_cw_repul_F0      = &
                colloid_get_cw_repul_F0(colloids,stat_info_sub)
           
           CALL colloid_get_shape(colloids, coll_shape,stat_info_sub)
           CALL colloid_get_radius(colloids, coll_radius,stat_info_sub)
           CALL colloid_get_freq(colloids, coll_freq,stat_info_sub)
           CALL colloid_get_m(colloids, coll_m,stat_info_sub)
           CALL colloid_get_mmi(colloids, coll_mmi,stat_info_sub)
           CALL colloid_get_x(colloids,coll_x,stat_info_sub)
           CALL colloid_get_v(colloids,coll_v,stat_info_sub)
           CALL colloid_get_accumulation_vector(colloids,&
                coll_acc_vector,stat_info_sub)
           CALL colloid_get_theta(colloids,coll_theta,stat_info_sub)
           CALL colloid_get_omega(colloids,coll_omega,stat_info_sub)
           
           WRITE(cbuf, '(A,I3)') 'coll_rho_type = ', coll_rho_type
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           WRITE(cbuf, '(A,L)') 'coll_translate = ', coll_translate
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           WRITE(cbuf, '(A,L)') 'coll_rotate = ', coll_rotate
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           WRITE(cbuf, '(A,I3)') 'coll_place = ', coll_place
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           WRITE(cbuf, '(A,I3)') 'coll_noslip = ', coll_noslip
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
              
           WRITE(cbuf, '(A,I3)') 'coll_body_force_type = ', coll_body_force_type
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           IF(num_dim == 2) THEN
              WRITE(cbuf, '(2(A,E16.8))') &
                   'coll_body_force = ',coll_body_force(1), ',', coll_body_force(2)
           ELSE IF (num_dim == 3 ) THEN
              WRITE(cbuf, '(3(A,E16.8))') &
                   'coll_body_force = ',coll_body_force(1), ',', coll_body_force(2), &
                   ',', coll_body_force(3)
           END IF
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '              
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
   
           WRITE(cbuf, '(A,E16.8)') 'cc_lub_cut_off = ', coll_cc_lub_cut_off
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           WRITE(cbuf, '(A,E16.8)') 'cc_lub_cut_on = ', coll_cc_lub_cut_on
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           WRITE(cbuf, '(A,E16.8)') 'cc_repul_cut_off = ', coll_cc_repul_cut_off
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
       
           WRITE(cbuf, '(A,E16.8)') 'cc_repul_cut_on = ', coll_cc_repul_cut_on
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
       
           WRITE(cbuf, '(A,E16.8)') 'cc_repul_F0 = ', coll_cc_repul_F0
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
       
           WRITE(cbuf, '(A,E16.8)') 'cw_lub_cut_off = ', coll_cw_lub_cut_off
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           WRITE(cbuf, '(A,E16.8)') 'cw_lub_cut_on = ', coll_cw_lub_cut_on
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
    
           WRITE(cbuf, '(A,E16.8)') 'cw_repul_cut_off = ', coll_cw_repul_cut_off
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
       
           WRITE(cbuf, '(A,E16.8)') 'cw_repul_cut_on = ', coll_cw_repul_cut_on
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           
           WRITE(cbuf, '(A,E16.8)') 'cw_repul_F0 = ', coll_cw_repul_F0
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
           WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '

           DO i = 1, num_colloid
              
              WRITE(cbuf, '(A,I3)') 'coll_shape = ', coll_shape(i)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
              
              IF (num_dim == 2) THEN
                 WRITE(cbuf, '(2(A,E16.8))') 'coll_radius = ', &
                      coll_radius(1,i), ',', coll_radius(2,i)
              ELSE
                 WRITE(cbuf, '(3(A,E16.8))') 'coll_radius = ', &
                      coll_radius(1,i), ',', coll_radius(2,i),&
                      ',', coll_radius(3,i)
              END IF
              
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
              
              
              WRITE(cbuf, '(A,I3)') 'coll_freq = ', coll_freq(i)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
              
              WRITE(cbuf, '(A,E16.8)') 'coll_m = ', coll_m(i)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
              
              WRITE(cbuf, '(3(A,E16.8))') &
                   'coll_mmi = ', coll_mmi(1,i), ',', coll_mmi(2,i), ',', coll_mmi(3,i)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
              
              IF(num_dim == 2) THEN
                 WRITE(cbuf, '(2(A,E16.8))') &
                      'coll_x = ',coll_x(1,i),',', coll_x(2,i)
              ELSE IF (num_dim == 3 ) THEN
                 WRITE(cbuf, '(3(A,E16.8))') &
                      'coll_x = ',coll_x(1,i),',', coll_x(2,i), ',', coll_x(3,i)
              END IF
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
              
              IF(num_dim == 2) THEN
                 WRITE(cbuf, '(2(A,E16.8))') &
                      'coll_v = ',coll_v(1,i),',', coll_v(2,i)
              ELSE IF (num_dim == 3 ) THEN
                 WRITE(cbuf, '(3(A,E16.8))') &
                      'coll_v = ',coll_v(1,i),',', coll_v(2,i), ',', coll_v(3,i)
              END IF
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
              
              WRITE(cbuf, '(4(A,E16.8))') &
                   'coll_rot_vector = ',coll_acc_vector(1,i),',', &
                   coll_acc_vector(2,i), ',', coll_acc_vector(3,i),',', coll_acc_vector(4,i)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
              
              
              WRITE(cbuf, '(3(A,E16.8))') 'coll_theta = ', coll_theta(1,i), &
                   ',', coll_theta(2,i), ',', coll_theta(3,i)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
              
              WRITE(cbuf, '(3(A,E16.8))') &
                   'coll_omega = ',coll_omega(1,i),',', coll_omega(2,i), ',', coll_omega(3,i)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
              WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
              
           END DO ! i = 1, num_colloid
           
        END IF ! num_specis > 1 && num_colloid > 0
        
        !----------------------------------------------------
        ! Write Boundary conditions of computational domain.
        !----------------------------------------------------
        
        CALL physics_get_bcdef(d_physics,bcdef,stat_info_sub)
        
        IF(num_dim == 2) THEN
           WRITE(cbuf, '(4(A,I2))') &
                'bcdef = ',bcdef(1),&
                (',', bcdef(i), i=2,4)
        ELSE IF (num_dim == 3 ) THEN
           WRITE(cbuf, '(6(A,I2))') &
                'bcdef = ',bcdef(1),&
                (',', bcdef(i), i=2,6)           
        END IF
        
        WRITE(UNIT=this%restart_physics_unit,&
             FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,&
             FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        IF(num_dim == 2) THEN
           WRITE(cbuf, '(4(A,I2))') &
                'shear_type = ',shear_type(1),&
                (',', shear_type(i),i=2,4)
        ELSE IF (num_dim == 3 ) THEN
           WRITE(cbuf, '(6(A,I2))') &
                'shear_type = ',shear_type(1),&
                (',', shear_type(i),i=2,6)
        END IF
        
        WRITE(UNIT=this%restart_physics_unit,&
             FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,&
             FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        
        IF( num_dim == 2 ) THEN
           WRITE(cbuf, '(4(A,E16.8))') &
                'shear_v = ',shear_v(2,1), &
                ',', shear_v(2,2), &
                ',', shear_v(1,3), &
                ',', shear_v(1,4)
           
           
        ELSE IF ( num_dim == 3 ) THEN           
           WRITE(cbuf, '(12(A,E16.8))') &
                'shear_v = ',shear_v(2,1),&
                ',',shear_v(3,1), &
                ',',shear_v(2,2),',',shear_v(3,2), &
                ',',shear_v(1,3),',',shear_v(3,3), &
                ',',shear_v(1,4),',',shear_v(3,4), &
                ',',shear_v(1,5),',',shear_v(2,5), &
                ',',shear_v(1,6),',',shear_v(2,6)
                 
        END IF
        
        WRITE(UNIT=this%restart_physics_unit,&
             FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,&
             FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        IF(num_dim == 2) THEN
           WRITE(cbuf, '(4(A,E16.8))') &
                'shear_freq = ',shear_freq(1),&
                (',', shear_freq(i),i=2,4 ) 
           
        ELSE IF (num_dim == 3 ) THEN
           WRITE(cbuf, '(6(A,E16.8))') &
                'shear_freq = ',shear_freq(1),&
                (',', shear_freq(i),i=2,6 ) 
        END IF
        
        
        WRITE(UNIT=this%restart_physics_unit,&
             FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,&
             FMT='(A)',IOSTAT=stat_info_sub) ' '
        
        WRITE(cbuf, '(A,I3)') 'wall_rho_type = ', wall_rho_type
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub) ' '
   
        
        WRITE(cbuf, '(A,I2)') 'wall_noslip = ', wall_noslip
        
        
        WRITE(UNIT=this%restart_physics_unit,&
             FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(UNIT=this%restart_physics_unit,&
             FMT='(A)',IOSTAT=stat_info_sub) ' '
        WRITE(UNIT=this%restart_physics_unit,&
             FMT='(A)',IOSTAT=stat_info_sub) ' '

        !------------------------------------------
        ! Write the ending of restart physics file.
        !------------------------------------------
        
        WRITE(cbuf, '(A)') '#---------------------------------------------'
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(cbuf, '(A)') '# The restart physics file ends here '
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        WRITE(cbuf, '(A)') '#---------------------------------------------'
        WRITE(UNIT=this%restart_physics_unit,FMT='(A)',IOSTAT=stat_info_sub)  TRIM(cbuf)
        
        
        IF( stat_info_sub /= 0 ) THEN
              
           PRINT *,"io_write_restart_physics : ",&
                "Writting into statis file failed!"
           
           stat_info = -1
           GOTO 9999
        END IF

        WRITE(cbuf,'(2A)') 'Restart physics config written to ',&
             TRIM(file_name)
        PRINT *, "***", TRIM(cbuf)
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(min_phys)) THEN
           DEALLOCATE(min_phys)
        END IF
        
        IF(ASSOCIATED(max_phys)) THEN
           DEALLOCATE(max_phys)
        END IF
        
        IF(ASSOCIATED(num_part_dim)) THEN
           DEALLOCATE(num_part_dim)
        END IF
        
        IF(ASSOCIATED(body_force)) THEN
           DEALLOCATE(body_force)
        END IF

        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF
        
        IF(ASSOCIATED(shear_type)) THEN
           DEALLOCATE(shear_type)
        END IF
        
        IF(ASSOCIATED(shear_v)) THEN
           DEALLOCATE(shear_v)
        END IF
        
        IF(ASSOCIATED(shear_freq)) THEN
           DEALLOCATE(shear_freq)
        END IF
        
        IF(ASSOCIATED(coll_body_force)) THEN
           DEALLOCATE(coll_body_force)
        END IF
        
        IF(ASSOCIATED(coll_shape)) THEN
           DEALLOCATE(coll_shape)
        END IF
        
        IF(ASSOCIATED(coll_radius)) THEN
           DEALLOCATE(coll_radius)
        END IF
        
        IF(ASSOCIATED(coll_m)) THEN
           DEALLOCATE(coll_m)
        END IF
        
        IF(ASSOCIATED(coll_mmi)) THEN
           DEALLOCATE(coll_mmi)
        END IF
        
        IF(ASSOCIATED(coll_x)) THEN
           DEALLOCATE(coll_x)
        END IF
        
        IF(ASSOCIATED(coll_v)) THEN
           DEALLOCATE(coll_v)
        END IF
        
        IF(ASSOCIATED(coll_acc_vector)) THEN
           DEALLOCATE(coll_acc_vector)
        END IF
        
        IF(ASSOCIATED(coll_theta)) THEN
           DEALLOCATE(coll_theta)
        END IF
        
        IF(ASSOCIATED(coll_omega)) THEN
           DEALLOCATE(coll_omega)
        END IF
        
        CLOSE(this%restart_physics_unit)
        
	RETURN 
 
      END SUBROUTINE io_write_restart_physics
