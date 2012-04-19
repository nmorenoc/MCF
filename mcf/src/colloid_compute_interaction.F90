      SUBROUTINE colloid_compute_interaction(this,&
           comm,MPI_PREC,drag,torque,FB,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_interaction
        !----------------------------------------------------
        !
        ! Purpose     : Compute colloid-colloid and
        !               colloid-wall interactions,
        !               if there is wall.
        !               The interactions may include:
        !               1)lubrication force correction,
        !               2)repulsive force to prevent
        !               overlap(contact force).
        !               One can also use 2) as physical
        !               DLVO repulsive force.
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     :
        !
        ! Revisions   : V0.1 19.11.2010, original version.
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
        ! Arguments
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(OUT)              :: this
        INTEGER, INTENT(IN)                     :: comm
        INTEGER, INTENT(IN)                     :: MPI_PREC
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: drag
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: torque
        REAL(MK), DIMENSION(:,:), INTENT(OUT)   :: FB
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim,num, dim2
        INTEGER                                 :: i,j        
        REAL(MK)                                :: cut_off
        REAL(MK)                                :: ghost_size
        REAL(MK), DIMENSION(3)                  :: min_phys_i
        REAL(MK), DIMENSION(3)                  :: max_phys_i
        REAL(MK), DIMENSION(:), POINTER         :: min_sub
        REAL(MK), DIMENSION(:), POINTER         :: max_sub
        REAL(MK), DIMENSION(3)                  :: min_sub_o
        REAL(MK), DIMENSION(3)                  :: max_sub_o
        
        LOGICAL                                 :: out        
        REAL(MK),DIMENSION(:,:),POINTER         :: x_t
        INTEGER, DIMENSION(:),  POINTER         :: sid_t
        INTEGER                                 :: num_t
        
        REAL(MK),DIMENSION(:,:),POINTER         :: x_ghost
        INTEGER, DIMENSION(:),  POINTER         :: sid_ghost
        INTEGER                                 :: num_g
        INTEGER                                 :: num_ghost
        
        LOGICAL                                 :: in
        REAL(MK),DIMENSION(:,:),POINTER         :: x_p
        INTEGER, DIMENSION(:), POINTER          :: sid_p
        REAL(MK),DIMENSION(3)                   :: F_ij
        REAL(MK),DIMENSION(3)                   :: F_i,F_j
        REAL(MK),DIMENSION(3)                   :: T_i,T_j
        REAL(MK),DIMENSION(:,:),POINTER         :: F, T
        REAL(MK),DIMENSION(:,:),POINTER         :: F_p, T_p
        REAL(MK),DIMENSION(:,:),POINTER         :: F_t, T_t
        INTEGER                                 :: num_p
        
        REAL(MK),DIMENSION(:,:), POINTER        :: x_p_ghost
        INTEGER, DIMENSION(:), POINTER          :: sid_p_ghost
        INTEGER                                 :: num_p_ghost

        REAL(MK), DIMENSION(3)                  :: length
        
        REAL(MK),DIMENSION(3,6)                 :: FB_lub
        REAL(MK),DIMENSION(3,6)                 :: FB_repul
        
        INTEGER                                 :: num_wall_solid
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        dim   = this%num_dim
        dim2  = dim * 2
        num   = this%num_colloid

        !----------------------------------------------------
        ! Set drag and torque of all colloids to zero
        ! before any operations.
        !----------------------------------------------------
        
        this%drag(1:dim,1:num) = 0.0_MK
        this%torque(1:3,1:num) = 0.0_MK
    
        IF ( SIZE(FB,1) /= dim ) THEN
           PRINT *, "colloid_compute_interaction : ",&
                "FB dimension do not match ! "
           stat_info = - 1
           GOTO 9999           
        END IF

        IF ( SIZE(FB,2) /= dim2 ) THEN
           PRINT *, "colloid_compute_interaction : ",&
                "FB number do not match ! "
           stat_info = - 1
           GOTO 9999           
        END IF
        
        FB(:,:) = 0
        
        NULLIFY(min_sub)
        NULLIFY(max_sub)
        
        NULLIFY(x_t)
        NULLIFY(sid_t)
        
        NULLIFY(x_ghost)
        NULLIFY(sid_ghost)
        
        NULLIFY(x_p)
        NULLIFY(sid_p)
        
        NULLIFY(F)        
        NULLIFY(F_p)
        NULLIFY(F_t)
        
        NULLIFY(T)        
        NULLIFY(T_p)
        NULLIFY(T_t)
      
        NULLIFY(x_p_ghost)
        NULLIFY(sid_ghost)
        
        length(1:dim) = &
             this%max_phys(1:dim) - this%min_phys(1:dim)
        
        num_wall_solid = &
             boundary_get_num_wall_solid(this%boundary,stat_info_sub)
        
        !----------------------------------------------------
        ! A : consider colloid-colloid lubrication correction
        !     contribution of total drag on colloids.
        ! B : consider colloid-colloid repulsive force
        !     (kind of contact force, or coating layer)
        !     contribution of total drag on colloids.
        ! 
        !     They are done in the same way.
        
        ! C : consider wall-colloid lubrication correction.
        ! D : consider wall-colloid repulsive force.
        !
        !     They are done in the same way.
        !----------------------------------------------------
        
        IF ( this%cc_lub_type == mcf_cc_lub_type_first .OR. &
             this%cc_repul_type >= mcf_cc_repul_type_Hookean .OR. &
             this%cw_lub_type == mcf_cw_lub_type_first .OR. &
             this%cw_repul_type >= mcf_cw_repul_type_Hookean ) THEN
           
           !-------------------------------------------------
           ! If we need any operation on colloids,
           ! we have to find which of them are in local
           ! sub domain.
           ! Get the boundary of this sub-domain.
           !-------------------------------------------------
           
           CALL technique_get_min_sub(this%tech,min_sub,stat_info_sub)
           CALL technique_get_max_sub(this%tech,max_sub,stat_info_sub)
      
           !-------------------------------------------------
           ! loop over each colloid in physics domain and
           ! record the ones, which are in this subdomain.
           !-------------------------------------------------
           
           ALLOCATE(x_p(dim,num))
           ALLOCATE(sid_p(num))
           num_p = 0
           
           DO j = 1, num
              
              !----------------------------------------------
              ! Assuming every one is inside the sub-domain
              ! first, but if one dimension is outside,
              ! it means it is outside.
              !----------------------------------------------

              in = .TRUE.
              
              DO i = 1, dim
                 
                 IF ( this%x(i,j) <  min_sub(i) .OR. &
                      this%x(i,j) >= max_sub(i) )  THEN
                    
                    in = .FALSE.
                    EXIT
                    
                 END IF
                 
              END DO
              
              !----------------------------------------------
              ! Record the one in this sub domain.
              !----------------------------------------------
              
              IF ( in ) THEN
                 
                 num_p            = num_p + 1
                 x_p(1:dim,num_p) = this%x(1:dim,j)
                 sid_p(num_p)     = j
                 
              END IF
              
           END DO ! j = 1, num
           
           !-------------------------------------------------
           ! Allocate force/torque memory of local real colloids.
           !-------------------------------------------------
           
           ALLOCATE(F_p(dim,num_p))
           F_p(1:dim,1:num_p) = 0.0_MK

           ALLOCATE(T_p(3,num_p))
           T_p(1:3,1:num_p) = 0.0_MK

           !-------------------------------------------------
           ! 1: If lubrication correction is required to amend
           !    forces between colloids, pair-wise forces
           !    are introduced for gap smaller than 
           !    cc_lub_cut_off.
           !
           ! 2D, 3D lurication theory formulations are different.
           !
           ! 2: if repusive force is needed,
           !    add it up to prevent overlaps between colloids.   
           !-------------------------------------------------
           
           IF ( this%cc_lub_type == mcf_cc_lub_type_first .OR. &
                this%cc_repul_type >= mcf_cc_repul_type_Hookean ) THEN
              
              !----------------------------------------------
              ! Set cut off, thereafter ghost zone size,
              ! which is used for building ghost zone.
              !----------------------------------------------
              
              cut_off = 0.0_MK
              
              IF ( this%cc_lub_type == mcf_cc_lub_type_first ) THEN
                 
                 cut_off = this%cc_lub_cut_off
                 
              END IF
              
              IF ( this%cc_repul_type >= mcf_cc_repul_type_Hookean .AND. &
                   this%cc_repul_cut_off > cut_off) THEN
                 
                 cut_off = this%cc_repul_cut_off
                 
              END IF
              
              ghost_size = 2.0_MK*this%radius(1,1) + cut_off
              
              !----------------------------------------------
              ! Set up the in-line box of physical domain
              ! to search for colloids.
              ! The ones outside the in-line box may have 
              ! images as ghosts of certain sub-domains.
              !
              ! Note that currently only periodic boundary
              ! is considered, not Lees-Edwards yet.
              !----------------------------------------------
              
              min_phys_i(1:dim) = this%min_phys(1:dim)
              max_phys_i(1:dim) = this%max_phys(1:dim)
              
              DO i = 1, dim
                 
                 IF ( this%bcdef(2*i-1) == ppm_param_bcdef_periodic) THEN
                    
                    min_phys_i(i) = min_phys_i(i) + ghost_size
                    
                 END IF
                 
                 IF ( this%bcdef(2*i) == ppm_param_bcdef_periodic) THEN
                    
                    max_phys_i(i) = max_phys_i(i) - ghost_size
                    
                 END IF
                 
              END DO ! i = 1, dim
              
              !----------------------------------------------
              ! loop all colloids and find the ones outside
              ! in-line box, of which images may be ghost
              ! colloids of some sum-domains.
              ! There are at most "num" of them.
              !----------------------------------------------
              
              ALLOCATE(x_t(dim,num))
              ALLOCATE(sid_t(num))
              num_t = 0
              
              DO j = 1, num
                 
                 !-------------------------------------------
                 ! Find the ones outside the in-line.
                 ! Assuming it is not outside, and search
                 ! each dimension coordinate, once found
                 ! one dimension outside, it is outside.
                 !-------------------------------------------
                 
                 out = .FALSE.
                 
                 DO i = 1, dim
                    
                    IF ( this%x(i,j) < min_phys_i(i) .OR. &
                         this%x(i,j) >= max_phys_i(i) )  THEN
                       
                       out = .TRUE.
                       EXIT
                       
                    END IF
                    
                 END DO ! i = 1, dim
              
                 !-------------------------------------------
                 ! Record the one in this out-line zone.
                 !-------------------------------------------
                 
                 IF ( out ) THEN
                    
                    num_t            = num_t + 1
                    x_t(1:dim,num_t) = this%x(1:dim,j)
                    sid_t(num_t)     = j
                    
                 END IF
                 
              END DO ! j = 1, num
              
              !----------------------------------------------
              ! Create ghost colloids using x_t(:,:)
              ! according to boundary conditions.
              ! For 2D,
              ! there are at most, num_t*2*2 of them.
              ! For 3D,
              ! there are at most, num_t*2*2*2 of them.
              !----------------------------------------------
           
              ALLOCATE(x_ghost(1:dim,num_t*2**dim))
              ALLOCATE(sid_ghost(num_t*2**dim))
              num_ghost = 0
              
              DO i = 1, dim
                 
                 num_g = num_ghost
                 
                 DO j = 1, num_ghost
                    
                    IF ( x_ghost(i,j) < min_phys_i(i) ) THEN
                       
                       num_g                = num_g + 1
                       x_ghost(1:dim,num_g) = x_ghost(1:dim,j)
                       x_ghost(i,num_g)     = x_ghost(i,num_g) + length(i)
                       sid_ghost(num_g)     = sid_ghost(j)
                       
                    ELSE IF ( x_ghost(i,j) >= max_phys_i(i) ) THEN
                       
                       num_g                = num_g + 1
                       x_ghost(1:dim,num_g) = x_ghost(1:dim,j)
                       x_ghost(i,num_g)     = x_ghost(i,num_g) - length(i)
                       sid_ghost(num_g)     = sid_ghost(j)
                       
                    END IF
                    
                 END DO ! j = 1, num_ghost
                 
                 
                 DO j = 1, num_t
                    
                    IF ( x_t(i,j) < min_phys_i(i)  ) THEN
                       
                       num_g                = num_g + 1
                       x_ghost(1:dim,num_g) = x_t(1:dim,j)
                       x_ghost(i,num_g)     = x_ghost(i,num_g) + length(i)
                       sid_ghost(num_g)     = sid_t(j)
                       
                    ELSE IF ( x_t(i,j) >= max_phys_i(i)  ) THEN
                       
                       num_g                = num_g + 1
                       x_ghost(1:dim,num_g) = x_t(1:dim,j)                    
                       x_ghost(i,num_g)     = x_ghost(i,num_g) - length(i)
                       sid_ghost(num_g)     = sid_t(j)
                       
                    END IF
                    
                 END DO ! j = 1, num_t
                 
                 num_ghost = num_g
                 
              END DO ! i = 1, dim
              
              !----------------------------------------------
              ! colloid-colloid interactions only on local
              ! process.
              !----------------------------------------------
              
              !----------------------------------------------
              ! Calculate its outside ghost zone boundary of
              ! this sub-domain.
              !----------------------------------------------
           
              min_sub_o(1:dim) = min_sub(1:dim) - ghost_size
              max_sub_o(1:dim) = max_sub(1:dim) + ghost_size
              
              ALLOCATE(x_p_ghost(dim,num+num_ghost))
              ALLOCATE(sid_p_ghost(num+num_ghost))
              num_p_ghost = 0
              
              !----------------------------------------------
              ! loop over each colloid in physics domain and
              ! record the ones, which are in this subdomain
              ! ghost zone.
              !----------------------------------------------
           
              DO j = 1, num
                 
                 !-------------------------------------------
                 ! Assuming every one is inside the subdomain+
                 ! ghost zone, but if one dimension is outside,
                 ! it means it is outside.
                 !-------------------------------------------
                 
                 in = .TRUE.
                 
                 DO i = 1, dim
                    
                    IF ( this%x(i,j) <  min_sub_o(i) .OR. &
                         this%x(i,j) >= max_sub_o(i) )  THEN
                       
                       in = .FALSE.
                       EXIT
                       
                    END IF
                    
                 END DO
              
                 !----------------------------------------------
                 ! Record the one in this sub-domain+ghost zone.
                 !----------------------------------------------
              
                 IF ( in ) THEN
                    
                    out = .FALSE.
                    
                    DO i = 1, dim
                       
                       IF ( this%x(i,j) <  min_sub(i) .OR. &
                            this%x(i,j) >= max_sub(i) )  THEN
                          
                          out = .TRUE.                          
                          EXIT
                          
                       END IF
                       
                    END DO
                    
                    !----------------------------------------
                    ! Record the one in ghost zone of
                    ! this sub domain.
                    !----------------------------------------
                    
                    IF ( out ) THEN
                    
                       num_p_ghost                  = num_p_ghost + 1
                       x_p_ghost(1:dim,num_p_ghost) = this%x(1:dim,j)
                       sid_p_ghost(num_p_ghost)     = j
                       
                    END IF
                    
                 END IF
                 
              END DO ! j = 1, num
           
              !----------------------------------------------
              ! Find ghost colloids outside of physics domain,
              ! (global ghosts colloids)
              ! which are inside the sub-domain's ghost zone.
              !----------------------------------------------
           
              DO j = 1, num_ghost
                 
                 !----------------------------------------------
                 ! Assuming every one is inside the subdomain+
                 ! ghost zone, but if one dimension is outside,
                 ! it means it is outside.
                 !----------------------------------------------
                 
                 in = .TRUE.
                 
                 DO i = 1, dim
                    
                    IF ( x_ghost(i,j) <  min_sub_o(i) .OR. &
                         x_ghost(i,j) >= max_sub_o(i) )  THEN
                       
                       in = .FALSE.
                       EXIT
                       
                    END IF
                    
                 END DO
                 
                 !----------------------------------------------
                 ! Record the one in this ghost zone.
                 ! (it can not be in the sub-domain anyway,
                 ! since it is outside physical domain).
                 !----------------------------------------------
                 
                 IF ( in ) THEN
                    
                    num_p_ghost                  = num_p_ghost + 1
                    x_p_ghost(1:dim,num_p_ghost) = x_ghost(1:dim,j)
                    sid_p_ghost(num_p_ghost)     = sid_ghost(j)
                    
                 END IF
                 
              END DO ! j = 1, num
              
                         
              !----------------------------------------------
              ! Loop each pair of real colloids.
              !----------------------------------------------
              
              DO i = 1, num_p-1
                 
                 DO j = i+1, num_p
                    
                    IF ( this%cc_lub_type == mcf_cc_lub_type_first ) THEN
                       
                       CALL colloid_compute_lubrication_cc(this,&
                            x_p(1:dim,i), x_p(1:dim,j),&
                            this%v(1:dim,sid_p(i),1),&
                            this%v(1:dim,sid_p(j),1),&
                            this%omega(1:3,sid_p(i),1),&
                            this%omega(1:3,sid_p(j),1),&
                            sid_p(i), sid_p(j),&
                            F_i(1:dim),F_j(1:dim),&
                            T_i(1:3),T_j(1:3),stat_info_sub)
                       
                       F_p(1:dim,i) = F_p(1:dim,i) + F_i(1:dim)
                       F_p(1:dim,j) = F_p(1:dim,j) + F_j(1:dim)
                       
                       T_p(1:3,i) = T_p(1:3,i) + T_i(1:3)
                       T_p(1:3,j) = T_p(1:3,j) + T_j(1:3)
                       
                    END IF
                    
                    IF ( this%cc_repul_type >= mcf_cc_repul_type_Hookean ) THEN
                       
                       CALL colloid_compute_repulsion_cc(this,&
                            x_p(1:dim,i),x_p(1:dim,j),&
                            sid_p(i),sid_p(j),F_ij(1:dim),stat_info_sub)
                       
                       F_p(1:dim,i) = F_p(1:dim,i) + F_ij(1:dim)
                       F_p(1:dim,j) = F_p(1:dim,j) - F_ij(1:dim)
                       
                    END IF
                    
                 END DO  ! j = i+1, num_p
                 
              END DO  ! i = 1, num_p -1
              
              !----------------------------------------------
              ! Loop each pair of real colloids in this 
              ! sub-domain and colloids in its ghost zone.
              !----------------------------------------------
              
              DO i = 1, num_p
                 
                 DO j = 1, num_p_ghost
                    
                    IF ( this%cc_lub_type == mcf_cc_lub_type_first ) THEN
                       
                       CALL colloid_compute_lubrication_cc(this,&
                            x_p(1:dim,i), x_p_ghost(1:dim,j),&
                            this%v(1:dim,sid_p(i),1), &
                            this%v(1:dim,sid_p_ghost(j),1),&
                            this%omega(1:3,sid_p(i),1),&
                            this%omega(1:3,sid_p_ghost(j),1),&
                            sid_p(i), sid_p_ghost(j), &
                            F_i(1:dim),F_j(1:dim),&
                            T_i(1:3),T_j(1:3),stat_info_sub)
                       
                       F_p(1:dim,i) = F_p(1:dim,i) + F_i(1:dim)
                       
                       T_p(1:3,i) = T_p(1:3,i) + T_i(1:3)
                       
                    END IF
                    
                    IF ( this%cc_repul_type >= mcf_cc_repul_type_Hookean ) THEN
                       
                       CALL colloid_compute_repulsion_cc(this,&
                            x_p(1:dim,i), x_p_ghost(1:dim,j),&
                            sid_p(i),sid_p_ghost(j),&
                            F_ij(1:dim),stat_info_sub)
                       
                       F_p(1:dim,i) = F_p(1:dim,i) + F_ij(1:dim)
                       
                    END IF
                    
                 END DO ! j = 1, num_p_ghost
                  
             END DO ! i = 1, num_p
              
           END IF ! cc_lub_type == 1 OR cc_repul_type == 1
           
           
           !-------------------------------------------------
           ! Consider colloid-wall interactions.
           !-------------------------------------------------
           
           IF ( num_wall_solid > 0 .AND. &
                ( this%cw_lub_type > mcf_cw_lub_type_no .OR. &
                this%cw_repul_type > mcf_cw_repul_type_no) ) THEN
              
              !----------------------------------------------
              ! Set total force from colloids on boundary
              ! to zero.
              !----------------------------------------------
              
              FB_lub(:,:)   = 0.0_MK
              FB_repul(:,:) = 0.0_MK
              
              !----------------------------------------------
              ! Loop over each colloid on this sub-domain.
              !----------------------------------------------
              
              DO i = 1, num_p
                 
                 !-------------------------------------------
                 ! If lubrication correction is required to 
                 ! amend forces between colloid and wall, 
                 ! forces are introduced for colloid-wall gap
                 ! smaller than cw_lub_cut_off.
                 !
                 ! 2D, 3D lurication theory formulations
                 ! are different.
                 !-------------------------------------------
                 
                 IF ( this%cw_lub_type > mcf_cw_lub_type_no ) THEN
                    
                    PRINT *, "colloid_compute_interaction: ", &
                         "calling lubrication_cw !"
                    
                    CALL colloid_compute_lubrication_cw(this,&
                         x_p(1:dim,i),this%v(1:dim,sid_p(i),1),sid_p(i),&
                         F_i(1:dim),FB(1:dim,1:dim2),stat_info_sub)
                    
                    IF ( stat_info_sub /= 0 ) THEN
                       PRINT *, "colloid_compute_interaction : ", &
                            "Computing lubrication cw failed !"
                       stat_info = -1
                       GOTO 9999
                    END IF
                    
                    F_p(1:dim,i) = F_p(1:dim,i) + F_i(1:dim)
                    
                    FB_lub(1:dim,1:dim2) = &
                         FB_lub(1:dim,1:dim2) + FB(1:dim,1:dim2)
                    
                 END IF ! cw_lub_type
                 
                 
                 !-------------------------------------------
                 ! Add up repulsive force to prevent overlaps
                 ! between wall and colloid.
                 !-------------------------------------------
                 
                 IF ( this%cw_repul_type > mcf_cw_repul_type_no ) THEN
                    
                    CALL colloid_compute_repulsion_cw(this,&
                         x_p(1:dim,i),sid_p(i),F_i(1:dim),&
                         FB(1:dim,1:dim2),stat_info_sub)
                    
                    IF ( stat_info_sub /= 0 ) THEN
                       PRINT *, "colloid_compute_interaction : ", &
                            "Computing repulsion cw failed !"
                       stat_info = -1
                       GOTO 9999
                    END IF
                    
                    F_p(1:dim,i) = F_p(1:dim,i) + F_i(1:dim)

                    FB_repul(1:dim,1:dim2) = &
                         FB_repul(1:dim,1:dim2) + FB(1:dim,1:dim2)
                    
                 END IF ! cc_repul_type
                 
              END DO  ! i = 1, num_p
              
              FB(1:dim,1:dim2) = &
                   FB_lub(1:dim,1:dim2) + FB_repul(1:dim,1:dim2)
              
           END IF ! cw_lub OR cw_repul
           

           !-------------------------------------------------
           ! Map the forces/torque on the correct colloids.
           ! sid_p(:) are indices for this procedure.
           !-------------------------------------------------
           
           ALLOCATE(F(dim,num))
           F(1:dim,1:num) = 0.0_MK
           ALLOCATE(F_t(dim,num))
           F_t(1:dim,1:num) = 0.0_MK
           
           DO i = 1, num_p
              F_t(1:dim, sid_p(i)) = F_p(1:dim,i)
           END DO
           
           ALLOCATE(T(3,num))
           T(1:3,1:num) = 0.0_MK
           ALLOCATE(T_t(3,num))
           T_t(1:3,1:num) = 0.0_MK
        
           DO i = 1, num_p
              T_t(1:3, sid_p(i)) = T_p(1:3,i)
           END DO
           
#ifdef __MPI
           
           !-------------------------------------------------
           ! In MPI context, collect all contributions 
           ! from each process.
           !-------------------------------------------------
           
           CALL MPI_ALLREDUCE (F_t(:,:),F(:,:), &
                SIZE(F),MPI_PREC,MPI_SUM,comm,stat_info_sub)
           
           !-------------------------------------------------
           ! Put colloid-colloid interactions on the total
           ! drag globally.
           !-------------------------------------------------
           
           this%drag(1:dim,1:num) = &
                drag(1:dim, 1:num) + F(1:dim,1:num)

           !-------------------------------------------------
           ! In MPI context, collect all contributions 
           ! from each process.
           !-------------------------------------------------
           
           CALL MPI_ALLREDUCE (T_t(:,:),T(:,:), &
                SIZE(T),MPI_PREC,MPI_SUM,comm,stat_info_sub)
           
           !-------------------------------------------------
           ! Put colloid-colloid interactions on the total
           ! torque globally.
           !-------------------------------------------------
           
           this%torque(1:3,1:num) = &
                torque(1:3, 1:num) + T(1:3,1:num)
#else
           
           !-------------------------------------------------
           ! Put colloid-colloid interactions on the total
           ! drag globally.
           !-------------------------------------------------
           
           this%drag(1:dim,1:num) = &
                drag(1:dim, 1:num) + F_t(1:dim,1:num)
           
           this%torque(1:3,1:num) = &
                torque(1:3, 1:num) + T_t(1:3,1:num)
        
#endif
           
        END IF ! cc_lub OR cc_repul OR cw_lub OR cw_repul
        
        
9999    CONTINUE
        
        
        IF(ASSOCIATED(min_sub)) THEN
           DEALLOCATE(min_sub)
        END IF
        
        IF(ASSOCIATED(max_sub)) THEN
           DEALLOCATE(max_sub)
        END IF
        
        IF(ASSOCIATED(x_t)) THEN
           DEALLOCATE(x_t)
        END IF
        
        IF(ASSOCIATED(sid_t)) THEN
           DEALLOCATE(sid_t)
        END IF
        
        IF(ASSOCIATED(x_ghost)) THEN
           DEALLOCATE(x_ghost)
        END IF
        
        IF(ASSOCIATED(sid_ghost)) THEN
           DEALLOCATE(sid_ghost)
        END IF
                
        IF(ASSOCIATED(x_p)) THEN
           DEALLOCATE(x_p)
        END IF
        
        IF(ASSOCIATED(sid_p)) THEN
           DEALLOCATE(sid_p)
        END IF
        
        IF(ASSOCIATED(F)) THEN
           DEALLOCATE(F)
        END IF

        IF(ASSOCIATED(F_p)) THEN
           DEALLOCATE(F_p)
        END IF

        IF(ASSOCIATED(F_t)) THEN
           DEALLOCATE(F_t)
        END IF

        IF(ASSOCIATED(T)) THEN
           DEALLOCATE(T)
        END IF

        IF(ASSOCIATED(T_p)) THEN
           DEALLOCATE(T_p)
        END IF

        IF(ASSOCIATED(T_t)) THEN
           DEALLOCATE(T_t)
        END IF
        
        IF(ASSOCIATED(x_p_ghost)) THEN
           DEALLOCATE(x_p_ghost)
        END IF
        
        IF(ASSOCIATED(sid_p_ghost)) THEN
           DEALLOCATE(sid_p_ghost)
        END IF
        
        RETURN          
        
      END SUBROUTINE colloid_compute_interaction
      
      
      
