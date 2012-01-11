      SUBROUTINE technique_init(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : techinique_init
        !----------------------------------------------------
        !
        ! Purpose     : Default construtor of technique
        !               class.
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
        TYPE(Technique),INTENT(OUT)             :: this 
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !----------------------------------------------------
        ! Init values using basic ppm functionality
        !----------------------------------------------------
        
        stat_info = 0
        
        NULLIFY(this%min_phys_t)
        NULLIFY(this%max_phys_t)
        NULLIFY(this%bcdef)
        
        this%igroup = 1
        this%ngroup = 1
        this%comm   = 0
        this%rank   = 0
        this%num_proc = 1          
        this%name_proc = "Azumi-2-Quadcores" 
        
        this%MPI_PREC = 0
        this%name_MPI_PREC = "NULL"
        
        this%decomp = ppm_param_decomp_bisection
        
        this%assig = ppm_param_assign_internal
        
        this%ghost_size   = 0.0_MK
        
        this%topo_id   = 1
        
        NULLIFY(this%min_sub)
        NULLIFY(this%max_sub)
        NULLIFY(this%sub_cost)
        NULLIFY(this%sub2proc)
        
        this%num_sub_tot = 0
        
        NULLIFY(this%sub_list)
        
        this%num_sub     = 0
        
        NULLIFY(this%sub_bcdef)

        NULLIFY(this%num_cell_dim_sub)
        ! Cell list by default
        this%neighbor_list = 2
        
        NULLIFY(this%cell_list)
        NULLIFY(this%inp)
        NULLIFY(this%jnp)        
        this%nnp = 0
        

        NULLIFY(this%num_cell_dim_sub_c)
        ! Cell list by default
        this%neighbor_list_c = 2
        
        NULLIFY(this%cell_list_c)
        NULLIFY(this%inp_c)
        NULLIFY(this%jnp_c)        
        this%nnp_c = 0

        this%ppm_debug    = 0
        
        this%ppm_log_unit  = 99    
        
        
9999    CONTINUE
        
        
        RETURN

      END SUBROUTINE technique_init
      
        
      SUBROUTINE technique_init_parallelization(this,igroup,ngroup,stat_info)
        
        !----------------------------------------------------
        ! Arguments :
        ! 
        ! igroup    : index of subgroup/communicator.
        ! ngroup    : total number of groups/communicators.
        !----------------------------------------------------
        
        TYPE(Technique),INTENT(OUT)     :: this
        INTEGER, INTENT(IN)             :: igroup
        INTEGER, INTENT(IN)             :: ngroup
        INTEGER, INTENT(OUT)            :: stat_info


        !----------------------------------------------------
        ! Local variables :
        !
        ! nproc  : total number of processes in MPI_COMM_WORLD
        ! mproc_group : maximum number of procs in new sub-group
        ! nproc_group : real number of procs in new sub-group
        ! ranks       : ranks in this group
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: nproc
        INTEGER                         :: mproc_group
        INTEGER                         :: nproc_group
        INTEGER,DIMENSION(:), POINTER   :: ranks
        
        INTEGER                         :: name_len, i
        INTEGER                         :: orig_group
        INTEGER                         :: new_group
        
        !--------------------------------
        ! Initialization
        !--------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        NULLIFY(ranks)
        
        !------------------------------------------
        ! Check the igroup and ngroup's validatiy.
        !------------------------------------------
        
        IF (igroup <=0 .OR. &
             ngroup <=0 .OR. &
             igroup > ngroup ) THEN
           PRINT *, "technique_init_parallelization : ", &
                "igroup or ngroup is wrong"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%igroup = igroup
        this%ngroup = ngroup
        
#ifdef __MPI
        
        !--------------------------------
        ! Initialize MPI.
        !--------------------------------    
        
        CALL MPI_Init(stat_info_sub)
        
        !--------------------------------
        ! Get total number of processors 
        ! for communicator MPI_COMM_WORLD
        !--------------------------------
        
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,stat_info_sub)
        
        !--------------------------------
        ! If desird subgroups are bigger
        ! than 1, we calculate number
        ! of processes in each sub-group.
        !--------------------------------
        
        IF( ngroup > 1) THEN 
           
           mproc_group = CEILING(REAL(nproc,MK)/REAL(ngroup,MK))
           nproc_group = mproc_group
           
           !-----------------------------
           ! Last sub-groups has the rest
           ! of processes, which might be
           ! smaller than mproc_group.
           !-----------------------------
           
           IF ( igroup == ngroup ) THEN
              
              nproc_group = &
                   nproc-(ngroup-1)*mproc_group              
              
           END IF
           
           !-----------------------------
           ! Generate an array which 
           ! contains all ranks that
           ! will be in same sub-group
           ! as this process.
           !-----------------------------
           
           ALLOCATE(ranks(nproc_group))
           
           DO i =1, nproc_group
              ranks(i) = (igroup-1) * mproc_group + i - 1
           END DO
           
           !-----------------------------
           ! Get the original group, which
           ! was for communicator
           ! MPI_COMM_WORLD.
           !-----------------------------
           
           CALL MPI_COMM_GROUP(MPI_COMM_WORLD,&
                orig_group,stat_info_sub)
           
           !-----------------------------
           ! Divide tasks into distinct
           ! subgroups.
           !-----------------------------
           
           CALL MPI_GROUP_INCL(orig_group,nproc_group,&
                ranks(1:nproc_group),new_group, stat_info_sub)
           
           
           !-----------------------------
           ! Create new communicator handel
           ! for each sub-group.
           !-----------------------------
           
           CALL MPI_COMM_CREATE(MPI_COMM_WORLD,new_group,&
                this%comm, stat_info_sub)
           
           
        ELSE
           
           this%comm = MPI_COMM_WORLD
           
        END IF
        
        
        !--------------------------------
        ! Get total number of processes
        ! in this communicator;
        ! rank of this processes in 
        ! this communicator.
        !--------------------------------
        
        
        CALL MPI_COMM_RANK(this%comm,this%rank,stat_info_sub)
        
        CALL MPI_COMM_SIZE(this%comm,this%num_proc,stat_info_sub)
        
        !--------------------------------
        ! Get name of processor
        !--------------------------------
        
        CALL MPI_GET_PROCESSOR_NAME(this%name_proc,name_len,stat_info_sub)
        
        !--------------------------------
        ! Define MPI precision
        !--------------------------------
        
        IF ( MK == ppm_kind_single ) THEN
           
           this%MPI_PREC = MPI_REAL
           this%name_MPI_PREC = 'Single Precision'
           
        ELSE IF ( MK == ppm_kind_double ) THEN
           
           this%MPI_PREC = MPI_DOUBLE_PRECISION
           this%name_MPI_PREC = 'Double Precision'
           
        ELSE
           
           WRITE(*,'(A,I4)')'Unknown Precision : ',MK
           stat_info_sub = -1
           
        END IF
        
        IF ( stat_info_sub /= 0 ) THEN
           
           stat_info = stat_info_sub
           PRINT *, "MPI Precision Error!"
           GOTO 9999
           
        END IF
#else
        
        this%rank = 0
        this%num_proc = 1          
        this%name_proc = "Azumi-OneCore"        
        
#endif
        
9999    CONTINUE
        
        IF (ASSOCIATED(ranks)) THEN
           DEALLOCATE(ranks)
        END IF
        
        RETURN
        
        
      END SUBROUTINE technique_init_parallelization
      
      
      SUBROUTINE technique_init_ppm(this,d_num_dim,&
           d_min_phys_t,d_max_phys_t,d_ghost_size,d_bcdef,stat_info)
        
        USE ppm_module_user_util
        USE ppm_module_user_io
        USE ppm_module_topo 
        
        TYPE(Technique),INTENT(OUT)     :: this
        INTEGER, INTENT(IN)             :: d_num_dim
        REAL(MK), DIMENSION(:)          :: d_min_phys_t
        REAL(MK), DIMENSION(:)          :: d_max_phys_t
        REAL(MK), INTENT(IN)            :: d_ghost_size
        INTEGER, DIMENSION(:)           :: d_bcdef
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        !Local variables
        INTEGER                         :: stat_info_sub
        INTEGER 			:: tolexp
        INTEGER                         :: i,j
        INTEGER                         :: isub
        REAL(MK), DIMENSION(3)          :: sub_min
        REAL(MK), DIMENSION(3)          :: sub_max

       
        this%num_dim   = d_num_dim
        
        IF(ASSOCIATED(this%min_phys_t)) THEN
           DEALLOCATE(this%min_phys_t)
        END IF
        ALLOCATE(this%min_phys_t(d_num_dim))
        this%min_phys_t(1:d_num_dim) = &
             d_min_phys_t(1:d_num_dim)
        
        IF(ASSOCIATED(this%max_phys_t)) THEN
           DEALLOCATE(this%max_phys_t)
        END IF
        ALLOCATE(this%max_phys_t(d_num_dim))
        this%max_phys_t(1:d_num_dim) = &
             d_max_phys_t(1:d_num_dim)
        
        this%ghost_size = d_ghost_size
        
        IF(ASSOCIATED(this%bcdef)) THEN
           DEALLOCATE(this%bcdef)
        END IF
        ALLOCATE(this%bcdef(2*d_num_dim))
        this%bcdef(1:2*d_num_dim) = &
             d_bcdef(1:2*d_num_dim)
        
        
        tolexp = INT(LOG10(EPSILON(this%ghost_size))) +1
        
        mcf_machine_zero = 10.0_MK**tolexp
        
        CALL ppm_init(d_num_dim,MK,tolexp,&
             this%comm,this%ppm_debug,&
             stat_info_sub,this%ppm_log_unit)
        
        
        IF (stat_info_sub /= 0) THEN
           PRINT * ,'technique_init_ppm',&
                'Init ppm library failed!'
           stat_info = -1
           GOTO 9999
        END IF
        
        !-------------------------------------------------------------
        ! Set ppm log file unit 
        !-------------------------------------------------------------
        
        CALL ppm_io_set_unit(6,0,this%ppm_log_unit,stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN
           PRINT *, 'technique_init_ppm : ',&
                'Set ppm library IO unit failed!'
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Create topology, actually calling
        ! ppm_topo_mkgeom_s / ppm_topo_mkgeom_d.
        !----------------------------------------------------
        
        
        CALL ppm_mktopo(this%decomp,this%assig,&
             d_min_phys_t(1:d_num_dim),&
             d_max_phys_t(1:d_num_dim), &
             d_bcdef(1:2*d_num_dim),&
             this%ghost_size,this%topo_id,&
             this%min_sub,this%max_sub,&
             this%sub_cost,this%sub2proc,this%num_sub_tot,&
             this%sub_list,this%num_sub,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN 
           PRINT *, "technique_init_ppm : ", &
                'Failed to creat topolgoy!'
           stat_info = stat_info_sub
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Get each subdomain for each rank.
        !----------------------------------------------------
        
        ALLOCATE(this%proc2sub(0:this%num_sub_tot-1))
        
        DO i = 1, this%num_sub_tot
           this%proc2sub(this%sub2proc(i)) = i
        END DO
        
        !----------------------------------------------------
        ! Check each subdomain on local process,
        ! decide wether it is on the boundary.
        !----------------------------------------------------

        ALLOCATE(this%sub_bcdef(2*d_num_dim,this%num_sub))
        this%sub_bcdef(:,:) = 0
        
        DO j = 1, this%num_sub
           
           isub = this%sub_list(j)
           
           sub_min(1:d_num_dim) =  &
                this%min_sub(1:d_num_dim,isub)
           sub_max(1:d_num_dim) =  &
                this%max_sub(1:d_num_dim,isub)
           
           DO i = 1, d_num_dim
              
              IF ( ABS(sub_min(i) - d_min_phys_t(i)) < &
                   mcf_machine_zero ) THEN
                 
                 this%sub_bcdef(2*i-1,j) = 1
                 
              END IF
              
              IF ( ABS(sub_max(i) - d_max_phys_t(i)) < &
                   mcf_machine_zero ) THEN
                 
                 this%sub_bcdef(2*i,j) = 1
                 
              END IF
              
           END DO ! i = 1, num_dim
           
        END DO ! j = 1, num_sub
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE technique_init_ppm
      
      
      SUBROUTINE technique_display_parameters(this,stat_info)
        
        TYPE(Technique),INTENT(IN)      :: this
        INTEGER,INTENT(OUT)             :: stat_info
        
        INTEGER                         :: num_dim, j
        INTEGER                         :: num_sub_tot
        
        
        num_dim     = this%num_dim
        num_sub_tot = this%num_sub_tot
        
        stat_info = 0
   
        PRINT *, '------------------Start------------------'
        PRINT *, '     Technique  parameters'
        PRINT *, '-----------------------------------------'
        
        PRINT *, "igroup           : ", this%igroup
        PRINT *, "ngroup           : ", this%ngroup
        PRINT *, "comm             : ", this%comm
        PRINT *, "rank             : ", this%rank
        PRINT *, "num_proc         : ", this%num_proc
        PRINT *, "name_proc        : ", TRIM(this%name_proc)
        PRINT *, "MPI_precision    : ", TRIM(this%name_MPI_PREC)
        PRINT *, "decomp           : ", this%decomp
        PRINT *, "assig            : ", this%assig
        PRINT *, "ghost_size       : ", this%ghost_size
        PRINT *, "topo_id          : ", this%topo_id
        PRINT *, "min_sub          : "
        DO j = 1, num_sub_tot
           PRINT *, this%min_sub(1:num_dim,j)
        END DO
        PRINT *, "max_sub          : "
        DO j = 1, num_sub_tot
           PRINT *, this%max_sub(1:num_dim,j)
        END DO
        
        PRINT *, "sub2proc         :"
        DO j = 1, SIZE(this%sub2proc,1)
           PRINT *, j, this%sub2proc(j)
        END DO
        
        PRINT *, "proc2sub         :"
        DO j = 0, SIZE(this%proc2sub,1)-1
           PRINT *, j, this%proc2sub(j)
        END DO
        
        PRINT *, "num_sub_tot      : ", num_sub_tot
        PRINT *, "num_sub          : ", this%num_sub
        
        PRINT *, "sub_list         : ", this%sub_list(1:this%num_sub)
        
        PRINT *, "sub_bcdef        : ", &
             this%sub_bcdef(1:2*num_dim,1:this%num_sub)
        
        SELECT CASE (this%neighbor_list)
           
        CASE (1)
           PRINT *, "neighbor_list    : ", "Verlet List"
           
        CASE (2)
           PRINT *, "neighbor_list    : ", "Cell List"
           
        END SELECT
        
        PRINT *, "ppm_debug        : ", this%ppm_debug
        !PRINT *, "ppm_log_unit     : ", this%ppm_log_unit
        
        PRINT *, '-------------------End-------------------'
        
        RETURN
        
      END SUBROUTINE technique_display_parameters
