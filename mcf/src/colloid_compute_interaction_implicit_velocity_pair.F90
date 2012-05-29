      SUBROUTINE colloid_compute_interaction_implicit_velocity_pair(this,&
           comm,MPI_PREC,step,dt,drag,torque,FB,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_interaction_implicit_
        !               velocity_pair
        !----------------------------------------------------
        !
        ! Purpose     : 1: Update velocities from SPH forces
        !                  using explicit scheme.
        !               2: Update velocties from 
        !                  lubrication-correction force using 
        !                  implicit splitting scheme for 
        !                  pair-wise colloids.
        !               3: Update velocities from pair-wise 
        !                  repulsive forces using explicit 
        !                  scheme.
        !
        ! Routines    :
        !
        ! References  : Shardlow T. SIAM J. Sci. Comput. 2003.
        !               Litvinov S. et al. J. Comput. Phys. 2010.
        !
        ! Remarks     : The interactions may include:
        !               1) SPH forces
        !               2) lubrication correction force
        !               3) repulsive force to prevent
        !                  overlap(contact force).
        !               One can also use 3) 2nd type as physical
        !               DLVO repulsive force.
        !
        !               1) Consider only serial run for the moment.
        !               2) body force is not considered.
        !
        ! Revisions   : V0.3 29.05.2012, sweeping of implicit pairwise
        !               velocity is improved to be adaptive, i.e.,
        !               number of sweeps is always decreased or increased
        !               by 2 times to have velocity difference between
        !               numbers of sweeping below (L2 norm)tolerance. 
        !
        !               V0.2 18.05.2012, shift 2nd repulsive force
        !               down to have clear zero force at 5*cut_off.
        !
        !               V0.1 09.05.2012, original version.
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
        INTEGER, INTENT(IN)                     :: step
        REAL(MK), INTENT(IN)                    :: dt
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: drag
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: torque
        REAL(MK), DIMENSION(:,:), INTENT(OUT)   :: FB
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local physics variables:
        !
        ! dim   : dimension
        ! num   : number of total colloids
        !
        ! hn_l,hm_l : lubrication correction cut off and on.
        ! F0,F1     : lubrication correction parameters.
        ! hn_r,hm_r : repulsive force cut off and on.
        ! F0_repul  : repulsive force magnitude.
        ! num_wall_solid : number of solid wall boundary.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: v_backup
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: v_sph
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: v_update
        REAL(MK), DIMENSION(:), ALLOCATABLE     :: v_diff,v_normal
        
        INTEGER                                 :: dim, num, dim2
        INTEGER                                 :: ent
        
        REAL(MK)                                :: hn_r, hm_r
        REAL(MK)                                :: F0_repul
       
        REAL(MK),DIMENSION(3)                   :: Fi_repul
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: F_repul
        REAL(MK),DIMENSION(3,6)                 :: FB_repul
      
        INTEGER                                 :: num_wall_solid
        
        !----------------------------------------------------
        ! Local numerical variables:
        !
        ! num_sweep   : number of sweeps for pair-wise implicit
        !               splitting scheme.
        !----------------------------------------------------
        
        REAL(MK)                                :: sweep_tolerance
        LOGICAL                                 :: sweep_adaptive
        
        INTEGER                                 :: num_sweep
        INTEGER                                 :: num_sweep1
        INTEGER                                 :: num_sweep2
        REAL(MK)                                :: error_decrease
        REAL(MK)                                :: error_increase
        REAL(MK)                                :: error_normal
        LOGICAL                                 :: sweep_decrease
        LOGICAL                                 :: sweep_increase
        
        !----------------------------------------------------
        ! num_sub_step: number of sub-steps for explicit 
        !               scheme.
        ! ai,aj       : radius of colloid i,j.
        ! aa          : ai+aj.
        ! Aij         : lubrication correction coefficient for i,j.
        ! Bij         : repulsive force coefficient for i,j.
        !
        ! x_image     : position of closest image of colloid j to i.
        ! v_image     : velocity of closest image of colloid j to i.
        !----------------------------------------------------
        
        INTEGER                                 :: num_sub_step, sub_step
        REAL(MK)                                :: dt_sub
        
        REAL(MK)                                :: ai, aj, aa
        REAL(MK)                                :: Bij
        
        REAL(MK), DIMENSION(3)                  :: x_image, v_image 
        REAL(MK), DIMENSION(3)                  :: rij, eij
        REAL(MK)                                :: r, h
          
        REAL(MK)                                :: fa, dt_f
        INTEGER                                 :: i, j, k, m        
        
        !----------------------------------------------------
        ! coll_k   : total kinetic energy of colloids
        ! coll_mom : total momentum of colloids
        !----------------------------------------------------
        
        REAL(MK)                                :: coll_k
        REAL(MK), DIMENSION(3)                  :: coll_mom
        CHARACTER(LEN=MAX_CHAR)                 :: fbuf

        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------

        stat_info     = 0
        stat_info_sub = 0
        
        FB(:,:) = 0.0_MK

        !----------------------------------------------------
        ! set up physics parameters.
        !----------------------------------------------------

        dim   = this%num_dim
        num   = this%num_colloid
        dim2  = dim * 2
        ent   = dim * num
        
        hn_r  = this%cc_repul_cut_off
        hm_r  = this%cc_repul_cut_on
        F0_repul = this%cc_repul_F0
        
        ALLOCATE(F_repul(dim,num))
        F_repul(:,:) = 0.0_MK
  
        num_wall_solid = &
             boundary_get_num_wall_solid(this%boundary,stat_info_sub)
        
        !----------------------------------------------------
        ! Back up previous velocities.
        !----------------------------------------------------
        
        ALLOCATE(v_backup(dim,num))
        v_backup(1:dim,1:num) = this%v(1:dim,1:num,1)
        
        !----------------------------------------------------
        ! First update velocity from  SPH forces contribution.
        !
        ! Check if the force is too big for current time step.
        !----------------------------------------------------
        
        DO i = 1, num

           fa = SQRT(DOT_PRODUCT(drag(1:dim,i),drag(1:dim,i)))/&
                this%m(i)
           
           IF ( fa > mcf_machine_zero ) THEN
              
              dt_f = this%adapt_t_coef*SQRT(this%h/fa)
              
              IF ( dt_f < dt_sub ) THEN
                 PRINT *, "drag is too big for colloid",&
                      i, drag(1:dim,i), dt_f
              END IF
              
           END IF
           
           this%v(1:dim,i,1) =  this%v(1:dim,i,1) + &
                dt * drag(1:dim,i) / this%m(i)
           
        END DO ! i = 1, num
        
        !----------------------------------------------------
        ! Back up velocity after SPH forces contributions.
        !----------------------------------------------------
        
        ALLOCATE(v_sph(dim,num))
        v_sph(1:dim,1:num) = this%v(1:dim,1:num,1)
        
        !----------------------------------------------------
        ! Do the pairwise sweep num_sweep times to update
        ! velocities using only lubrcation correction forces.
        ! Time step at each sweep is dt_sweep = dt / num_sweep.
        !----------------------------------------------------
        
        IF ( this%cc_lub_type > mcf_cc_lub_type_no ) THEN
           
           !-------------------------------------------------
           ! set up parameters for implicit sweeps.
           !
           ! sweep_adaptive : check if sweeping is adaptive.
           ! num_sweep      : initial number of sweeps.
           !-------------------------------------------------
           
           sweep_adaptive  = &
                this%implicit_pair_sweep_adaptive
           num_sweep = this%implicit_pair_num_sweep
           
           !-------------------------------------------------
           ! Compute implicit veloicyt pairwise  interaction
           ! num_sweep times.
           !-------------------------------------------------
           
           CALL colloid_compute_interaction_implicit_velocity_pair_sweep(&
                this, dt, num_sweep, stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, __FILE__, __LINE__, &
                   "velocity pair sweeps failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Perform adapting number of sweeps, if sweeping
           ! is required adaptive.
           !-------------------------------------------------
           
           IF ( sweep_adaptive ) THEN
              
              !----------------------------------------------
              ! Get tolerance of between two numbers of sweeping.
              !----------------------------------------------
              
              sweep_tolerance = &
                   this%implicit_pair_sweep_tolerance

              !----------------------------------------------
              ! Back up velocity from previous num_sweep sweeps.
              !----------------------------------------------
              
              ALLOCATE(v_update(dim,num))
              v_update(1:dim,1:num) = this%v(1:dim,1:num,1)
              
              !----------------------------------------------
              ! Allocate memory for velocity difference.
              !----------------------------------------------

              ALLOCATE(v_diff(ent))
              ALLOCATE(v_normal(ent))
              
              !----------------------------------------------
              ! Set initial error to be an artibrarily 
              ! small number, and perform decreasing number
              ! of sweeps.
              !----------------------------------------------
              
              error_decrease = 0.0_MK
              num_sweep1     = num_sweep / 2
              sweep_decrease = .FALSE.
              
              !----------------------------------------------
              ! Minimum number of sweeps is 1.
              !----------------------------------------------

              DO  WHILE ( num_sweep1 >= 1 .AND. &
                   error_decrease < sweep_tolerance ) 
                 
                 !-------------------------------------------
                 ! Set velocity to v_sph, i.e., after
                 ! SPH forces contributions and perform
                 ! pairwise sweeping num_sweep1 times.
                 !-------------------------------------------
                 
                 this%v(1:dim,1:num,1) = v_sph(1:dim,1:num)
                 
                 CALL colloid_compute_interaction_implicit_velocity_pair_sweep(&
                      this, dt, num_sweep1, stat_info_sub)
                 
                 IF ( stat_info_sub /= 0 ) THEN
                    PRINT *, __FILE__, __LINE__, &
                         "velocity pair sweeps failed!"
                    stat_info = -1
                    GOTO 9999
                 END IF
                 
                 !-------------------------------------------
                 ! Calculate velocity different between
                 ! num_sweep and num_sweep1 sweeps.
                 !-------------------------------------------
                 
                 DO i = 1, dim
                    
                    v_diff((i-1)*num+1:i*num) = &
                         this%v(i,1:num,1) - v_update(i,1:num)
                    v_normal((i-1)*num+1:i*num) = &
                         v_update(i,1:num)
                    
                 END DO
                 
                 
                 !-------------------------------------------
                 ! Calculate error(L2 norm) beween two numbers 
                 ! of sweeping.
                 !-------------------------------------------
                 
                 error_decrease = &
                      tool_L2_norm(this%tool,v_diff(1:ent),stat_info_sub)
                 error_normal   = &
                      tool_L2_norm(this%tool,v_normal(1:ent),stat_info_sub)
                 
                 IF ( error_normal < mcf_machine_zero ) THEN
                    
                    error_normal = mcf_machine_zero
                    
                 END IF
                 
                 !PRINT *, "sweep1, error_normal: ", error_normal
                 
                 error_decrease = error_decrease / error_normal
                 
                 IF ( error_decrease < sweep_tolerance ) THEN
                    
                    !----------------------------------------
                    ! If error is small enough, decrease
                    ! number of sweeps again.
                    !----------------------------------------
                    
                    num_sweep   = num_sweep1
                    num_sweep1  = num_sweep / 2
                    v_update(1:dim,1:num) = this%v(1:dim,1:num,1)
                    sweep_decrease = .TRUE.
                    
                    !----------------------------------------
                    ! If number of sweeps was decreased,
                    ! save the current error.
                    !----------------------------------------
                    
                    this%implicit_pair_sweep_error = error_decrease
                    
                 END IF ! error_decrease < tolerance
                 
              END DO ! num_sweep1 > 0 AND error_decrease < tolerance
              
            
              IF ( sweep_decrease ) THEN
                 
                 !-------------------------------------------
                 ! If number of sweeps was decreased,
                 ! there is no need to perform increasing and
                 ! save the current sutiable number of sweeps.
                 !-------------------------------------------
                 
                 this%implicit_pair_num_sweep   = num_sweep     
                 
              ELSE
                 
                 !-------------------------------------------
                 ! If number of sweeps was not decreased,
                 ! we need to increase the number of sweeps.
                 ! There may be two reasons that we did not
                 ! decrease number of sweeps.
                 ! 1: num_sweep1 < 1 or
                 ! 2: error_decrease > sweep_tolerance.
                 !-------------------------------------------

                 sweep_increase = .FALSE.
                 num_sweep2 = num_sweep * 2
                 
                 IF ( num_sweep1 >= 1 ) THEN
                    
                    !----------------------------------------
                    ! set initial error_increase to be
                    ! error_decrease.
                    !----------------------------------------
                    
                    error_increase = error_decrease
                    
                 ELSE
                    
                    !----------------------------------------
                    ! set initial error to be an artibrarily 
                    ! big number.
                    !----------------------------------------
             
                    error_increase = 1.0e2_MK
                    
                 END IF
                 
                 !-------------------------------------------
                 ! perform number of sweeps until maximum 
                 ! number or error is small enough.
                 !-------------------------------------------

                 DO  WHILE ( num_sweep2 <= &
                      mcf_cc_lub_implicit_velocity_sweep_max .AND. &
                      error_increase > sweep_tolerance ) 
                    
                    !----------------------------------------
                    ! Set velocity to v_sph, i.e., after
                    ! SPH forces contributions and num_sweep2
                    ! times sweeping.
                    !----------------------------------------
                    
                    this%v(1:dim,1:num,1) = v_sph(1:dim,1:num)
                    
                    CALL colloid_compute_interaction_implicit_velocity_pair_sweep(&
                         this, dt, num_sweep2, stat_info_sub)
                    
                    IF ( stat_info_sub /= 0 ) THEN
                       PRINT *, __FILE__, __LINE__, &
                            "velocity pair sweeps failed!"
                       stat_info = -1
                       GOTO 9999
                    END IF
                    
                    !-------------------------------------------
                    ! Calculate velocity different between
                    ! num_sweep and num_sweep1 sweeps.
                    !-------------------------------------------
                    
                    DO i = 1, dim
                       
                       v_diff((i-1)*num+1:i*num) = &
                            this%v(i,1:num,1) - v_update(i,1:num)
                       v_normal((i-1)*num+1:i*num) = &
                            v_update(i,1:num)
                       
                    END DO
                    
                    !----------------------------------------
                    ! Calculate error beween two number of 
                    ! sweeps
                    !----------------------------------------
                    
                    error_increase = &
                            tool_L2_norm(this%tool,v_diff(1:ent),stat_info_sub)
                    
                    error_normal   = &
                         tool_L2_norm(this%tool,v_normal(1:ent),stat_info_sub)
                    
                    IF ( error_normal < mcf_machine_zero ) THEN
                       
                       error_normal = mcf_machine_zero
                       
                    END IF
                    
                    !PRINT *, "sweep2, error_normal: ", error_normal
                    
                    error_increase = error_increase / error_normal
                    
                    IF ( error_increase > sweep_tolerance ) THEN
                       
                       num_sweep  = num_sweep2
                       num_sweep2 = num_sweep * 2
                       v_update(1:dim,1:num) = this%v(1:dim,1:num,1)
                       sweep_increase = .TRUE.
                       
                    END IF ! error_increase > tolerance
                    
                 END DO ! num_sweep2 < sweep_max AND error_increase > tolerance
                 
                 !-------------------------------------------
                 ! No matter if it is increased, number of
                 ! sweeps and error have to be recorded.
                 !-------------------------------------------                 
          
                 this%implicit_pair_num_sweep   = num_sweep                 
                 this%implicit_pair_sweep_error = error_increase
                 
                 !-------------------------------------------
                 ! If not increase, velocity has to be restored.
                 !-------------------------------------------

                 IF (  .NOT. sweep_increase ) THEN
                    
                    this%v(1:dim,1:num,1) = v_update(1:dim,1:num)
                    
                 END IF
                 
              END IF ! NOT sweep_decrease

           END IF ! sweep_adaptive
           
           !PRINT *, step, this%implicit_pair_num_sweep, error
           
        END IF ! cc_lub_type > mcf_cc_lub_type_no
        
        !----------------------------------------------------
        ! Loop over each colloid and its interaction with 
        ! other colloids and walls,
        ! calculate the repulsive force between 
        ! colloid-colloid and colloid-wall.
        !
        ! Update velocity explicitly from  repulsive force
        ! between colloid-colloid and colloid-wall
        ! by two sub steps, each one dt_sub/2.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! set up parameters for explicit update.
        !----------------------------------------------------
        
        num_sub_step = this%explicit_sub_time_step
        dt_sub       = dt / num_sub_step
        
        !----------------------------------------------------
        ! Repulsive forces between colloid-colloid and
        ! colloid-wall.
        !----------------------------------------------------
        
        IF ( this%cc_repul_type > mcf_cc_repul_type_no .OR. &
             this%cw_repul_type > mcf_cw_repul_type_no ) THEN
           
           !-------------------------------------------------
           ! For now, we assume that each radius is the same.
           !-------------------------------------------------
           
           ai = this%radius(1,1)
           aj = ai
           aa = ai+aj
    
           !-------------------------------------------------
           ! As we use modified velocity verlet, dt_sub is
           ! splitted into two sub-steps and the force
           ! needs to be calculated only once per dt_sub.
           !
           ! Therefore, if sub_step > 1, first F_repul does 
           ! not need to be caulcated again, since the position
           ! is not changed.
           ! Only the second F_repul is calculated at 
           ! every sub_step.
           !-------------------------------------------------
           
           F_repul(1:dim,1:num) = 0.0_MK
           
           DO sub_step = 1, num_sub_step
              
              IF ( this%cc_repul_type > mcf_cc_repul_type_no .AND. &
                   sub_step == 1 ) THEN
                 
                 !-------------------------------------------
                 ! Repulsive force between colloid i,j.
                 !-------------------------------------------

                 DO i = 1, num - 1
                    
                    DO j = i + 1, num
                       
                       !-------------------------------------
                       ! Calculate gap and unit vector 
                       ! joining two colloids, images of 
                       ! colloids have to be considered 
                       ! according to different boundaries.
                       !-------------------------------------
                 
                       CALL colloid_nearest_image(this,&
                            this%x(1:dim,i),j, &
                            x_image(1:dim),rij(1:dim), &
                            v_image(1:dim),stat_info_sub)
                       
                       r  = SQRT(DOT_PRODUCT(rij(1:dim), rij(1:dim)))
                       h  = r - aa
                       eij(1:dim) = rij(1:dim) / r
                    
                       !-------------------------------------
                       ! If gap is smaller than 5*hn_r, 
                       ! it may need repulsive force.
                       !-------------------------------------
                       
                       IF ( h < 5.0_MK * hn_r ) THEN
                          
                          Bij = 0.0_MK
                          
                          SELECT CASE ( this%cc_repul_type )
                             
                          CASE ( mcf_cc_repul_type_Hookean )
                             
                             !-------------------------------
                             ! For linear spring force, it 
                             ! has clear zero value at hn_r.
                             !-------------------------------
                          
                             IF ( h < hn_r ) THEN
                                
                                !----------------------------
                                ! If gap smaller than minimal
                                ! allowed gap, set it to 
                                ! the pre-set minimum.
                                !----------------------------
                                
                                IF ( h < hm_r ) THEN
                                   
                                   h = hm_r
                                   
                                END IF
                                
                                Bij = F0_repul - F0_repul*h/hn_r
                                
                             END IF ! h < hn_r
                          
                          CASE ( mcf_cc_repul_type_DLVO )
                          
                             !-------------------------------
                             ! For DLVO force, it does not 
                             ! have clear zero at hn_r.
                             ! But at 5*hn_r, its value is 
                             ! smaller than F0/100.
                             ! Also minus the value at 5*hn_r
                             ! to have a clear zero value of
                             ! repulsive force at 5*hn_r
                             !-------------------------------
                             
                             !-------------------------------
                             ! If gap smaller than minimal 
                             ! allowed gap,
                             ! set it to the pre-set minimum.
                             !-------------------------------
                             
                             IF ( h < hm_r ) THEN
                                
                                h = hm_r
                                
                             END IF
                             
                             
                             Bij = F0_repul / hn_r * &
                                  EXP(-h/hn_r) / (1.0_MK-EXP(-h/hn_r)) - &
                                  F0_repul / hn_r * &
                                  EXP(-5.0_MK) / (1.0_MK-EXP(-5.0_MK))
                             
                          CASE DEFAULT
                             
                             PRINT *, __FILE__, __LINE__, &
                                  "no such repulsive force!"
                             stat_info = -1
                             GOTO 9999
                             
                          END SELECT ! cc_repul_type
                          
                          !----------------------------------
                          ! Save pair-wise repulsive forces
                          !----------------------------------
                          
                          F_repul(1:dim,i) = F_repul(1:dim,i) + &
                               Bij * eij(1:dim)
                          F_repul(1:dim,j) = F_repul(1:dim,j) - &
                               Bij * eij(1:dim)
                          
                       END IF ! h < hn_r
                       
                    END DO ! j = i +1, num
                    
                 END DO ! i = 1, num - 1
                 
              END IF ! cc_repul_type > mcf_cc_repul_type_no AND sub_step == 1
              
              !----------------------------------------------
              ! Repulsive forces between colloid-wall.
              !----------------------------------------------
              
              IF ( this%cw_repul_type > mcf_cw_repul_type_no .AND. &
                   sub_step == 1 ) THEN
                 
                 DO i = 1, num
                    
                    CALL colloid_compute_repulsion_cw(this,&
                         this%x(1:dim,i),i,Fi_repul(1:dim),&
                         FB_repul(1:dim,1:dim2),stat_info_sub)
                    
                    IF ( stat_info_sub /= 0 ) THEN
                       PRINT *, __FILE__, __LINE__, &
                            "calling repulsion cw failed !"
                       stat_info = -1
                       GOTO 9999
                    END IF
                    
                    F_repul(1:dim,i) = &
                         F_repul(1:dim,i) + Fi_repul(1:dim)
                    
                 END DO ! i = 1, num
                 
              END IF ! cw_repul_type > mcf_cw_repul_type_no AND sub_step == 1
              
              !----------------------------------------------
              ! Update velocity using repulsive force 
              ! by dt_sub/2.
              ! Check if the repulsive force is too big
              ! force dt_sub/2.
              !----------------------------------------------
              
              DO i = 1, num
                 
                 fa = SQRT(DOT_PRODUCT(F_repul(1:dim,i),F_repul(1:dim,i)))/&
                      this%m(i)
                 
                 IF ( fa > mcf_machine_zero ) THEN
                    
                    dt_f = this%adapt_t_coef*SQRT(this%h/fa)
                    
                    IF ( dt_f < dt_sub/2.0_MK ) THEN
                       
                       PRINT *, "First dt_sub: F_repul is too big for colloid",&
                            i, F_repul(1:dim,i), dt_f, sub_step
                       
                    END IF
                    
                 END IF
                 
                 this%v(1:dim,i,1) = &
                      this%v(1:dim,i,1) + &
                      0.5_MK *  F_repul(1:dim,i) * dt_sub / this%m(i)
                 
              END DO ! i = 1, num
              
              !----------------------------------------------
              ! Update position by dt_sub.
              !----------------------------------------------
              
              DO i = 1, num
                 
                 this%x(1:dim,i) = &
                      this%x(1:dim,i) + this%v(1:dim,i,1) * dt_sub
                 
              END DO ! i = 1, num
              
              !----------------------------------------------
              ! Adjust colloids position according to 
              ! boundary conditions.
              !----------------------------------------------
              
              CALL colloid_adjust_colloid(this,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 
                 PRINT *, __FILE__, __LINE__, &
                      "adjusting colloid failed!"
                 stat_info = -1
                 GOTO 9999
                 
              END IF
              
              !----------------------------------------------
              ! As some colloids may cross boundaries,
              ! compute theire new images(position and velocity).
              !----------------------------------------------
              
              CALL colloid_compute_image(this,stat_info_sub)
              
              IF ( stat_info_sub /=0 ) THEN
                 
                 PRINT *, __FILE__, __LINE__, &
                      "colloid computing image failed!"
                 stat_info = -1
                 GOTO 9999
                 
              END IF
              
              !----------------------------------------------
              ! Calculate the pairwise repulsive force again 
              ! using new positions.
              !----------------------------------------------
              
              F_repul(1:dim,1:num) = 0.0_MK
              
              IF ( this%cc_repul_type > mcf_cc_repul_type_no ) THEN
                 
                 DO i = 1, num - 1
                    
                    DO j = i + 1, num
                       
                       !-------------------------------------
                       ! Calculate gap and unit vector joining
                       ! two colloids, images of colloids have
                       ! to be considered according to different
                       ! boundary conditions.
                       !-------------------------------------
                       
                       CALL colloid_nearest_image(this,&
                            this%x(1:dim,i),j, &
                            x_image(1:dim),rij(1:dim), &
                            v_image(1:dim),stat_info_sub)
                       
                       r  = SQRT(DOT_PRODUCT(rij(1:dim), rij(1:dim)))
                       h  = r - aa
                       eij(1:dim) = rij(1:dim) / r
                       
                       !-------------------------------------
                       ! If gap is smaller than 5*hn_r, it
                       ! may need repulsive force.
                       !-------------------------------------
                       
                       IF ( h < 5.0_MK * hn_r ) THEN
                          
                          Bij = 0.0_MK
                          
                          SELECT CASE ( this%cc_repul_type )
                             
                          CASE ( mcf_cc_repul_type_Hookean )
                             
                             !-------------------------------
                             ! For linear spring force, it has 
                             ! clear zero at hn_r.
                             !-------------------------------
                             
                             IF ( h < hn_r ) THEN
                                
                                !----------------------------
                                ! If gap smaller than minimal 
                                ! allowed gap, 
                                ! set it to the pre-set minimum.
                                !----------------------------
                             
                                IF ( h < hm_r ) THEN
                                
                                   h = hm_r
                                   
                                END IF
                             
                                Bij = F0_repul - F0_repul*h/hn_r
                                
                             END IF ! h < hn_r
                             
                          CASE ( mcf_cc_repul_type_DLVO )
                             
                             !-------------------------------
                             ! For DLVO force, it does not have 
                             ! clear zero at hn. But at 5*hn, 
                             ! its value is smaller than F0/100.
                             ! Also minus the value at 5*hn_r
                             ! to have a clear zero value of
                             ! repulsive force at 5*hn_r                       
                             !-------------------------------
                          
                             !-------------------------------
                             ! If gap smaller than minimal 
                             ! allowed gap,
                             ! set it to the pre-set minimum.
                             !-------------------------------
                          
                             IF ( h < hm_r ) THEN
                                
                                h = hm_r
                                
                             END IF
                             
                             Bij = F0_repul / hn_r * &
                                  EXP(-h/hn_r) / (1.0_MK-EXP(-h/hn_r)) - &
                                  F0_repul / hn_r * &
                                  EXP(-5.0_MK) / (1.0_MK-EXP(-5.0_MK))
                             
                          CASE DEFAULT
                             
                             PRINT *, __FILE__, __LINE__, &
                                  "no such repulsive force!"
                             stat_info = -1
                             GOTO 9999
                             
                          END SELECT ! cc_repul_type
                          
                          !-----------------------------------
                          ! Save pair-wise repulsive forces
                          !----------------------------------
                          
                          F_repul(1:dim,i) = F_repul(1:dim,i) + &
                               Bij * eij(1:dim)
                          F_repul(1:dim,j) = F_repul(1:dim,j) - &
                               Bij * eij(1:dim)
                          
                       END IF ! h < hn_r
                       
                    END DO ! j = i + 1, num
                    
                 END DO ! i = 1, num - 1
                 
              END IF ! cc_repul_type > mcf_cc_repul_type_no
              
              !----------------------------------------------
              ! Repulsive force between colloid-wall.
              ! Save the repulsive force on the wall also.
              !----------------------------------------------
              
              FB(1:dim,1:dim2) = 0.0_MK
              
              IF ( this%cw_repul_type > mcf_cw_repul_type_no ) THEN
                 
                 DO i = 1, num
                    
                    CALL colloid_compute_repulsion_cw(this,&
                         this%x(1:dim,i),i,Fi_repul(1:dim),&
                         FB_repul(1:dim,1:dim2),stat_info_sub)
                    
                    IF ( stat_info_sub /= 0 ) THEN
                       PRINT *, __FILE__, __LINE__, &
                            "calling repulsion cw failed !"
                       stat_info = -1
                       GOTO 9999
                    END IF
                    
                    F_repul(1:dim,i) = &
                         F_repul(1:dim, i) + Fi_repul(1:dim)
                    FB(1:dim,1:dim2) = &
                         FB(1:dim,1:dim2) + FB_repul(1:dim,1:dim2)
                    
                 END DO ! i = 1, num
              
              END IF ! cw_repul_type > no
           
              !----------------------------------------------------
              ! Update velocity using repulsive force at 
              ! second dt_sub/2.
              ! Check if the repulsive force is too big
              ! force dt_sub/2.
              !----------------------------------------------------
              
              DO i = 1, num
                 
                 fa = SQRT(DOT_PRODUCT(F_repul(1:dim,i),F_repul(1:dim,i)))/&
                      this%m(i)
                 
                 IF ( fa > mcf_machine_zero ) THEN
                    
                    dt_f = this%adapt_t_coef*SQRT(this%h/fa)
                    
                    IF ( dt_f < dt_sub/2.0_MK ) THEN
                       
                       PRINT *, "Second dt_sub: F_repul is too big for colloid",&
                            i, F_repul(1:dim,i), dt_f,sub_step
                       
                    END IF
                    
                 END IF
                 
                 this%v(1:dim,i,1) = &
                      this%v(1:dim,i,1) + &
                      0.5_MK * F_repul(1:dim,i) * dt_sub /&
                      this%m(i)           
                 
              END DO ! i = 1, num
              
           END DO ! sub_step = 1, num_sub_step
           
        ELSE
           
           !-------------------------------------------------
           ! Without any repulsive force, update position by dt.
           !-------------------------------------------------
           
           DO i = 1, num
              
              this%x(1:dim,i) = &
                   this%x(1:dim,i) + this%v(1:dim,i,1) * dt
              
           END DO ! i = 1, num
           
           !-------------------------------------------------
           ! Adjust colloids position according to 
           ! boundary conditions.
           !-------------------------------------------------
              
           CALL colloid_adjust_colloid(this,stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              
              PRINT *, __FILE__, __LINE__, &
                   "adjusting colloid failed!"
              stat_info = -1
              GOTO 9999
              
           END IF
           
           !-------------------------------------------------
           ! As some colloids may cross boundaries,
           ! compute theire new images(position and velocity).
           !-------------------------------------------------
           
           CALL colloid_compute_image(this,stat_info_sub)
           
           IF ( stat_info_sub /=0 ) THEN
              
              PRINT *, __FILE__, __LINE__, &
                   "colloid computing image failed!"
              stat_info = -1
              GOTO 9999
              
           END IF
           
        END IF ! cc_repul_type > no AND cw_repul_type > no
        
        !----------------------------------------------------
        ! Work out the effective drag by calculating the
        ! difference between new and old velocity and
        ! assuming an Euler scheme was used as integrator.
        !----------------------------------------------------
        
        DO i = 1, num
           
           this%drag(1:dim,i) = &
                (this%v(1:dim,i,1) - v_backup(1:dim,i))*this%m(i) / dt
           
        END DO ! i = 1, num
        
        !----------------------------------------------------
        ! Update total torque on each colloid
        !----------------------------------------------------
        
        this%torque(1:3,1:num) = torque(1:3,1:num)
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE colloid_compute_interaction_implicit_velocity_pair
      
      
      
