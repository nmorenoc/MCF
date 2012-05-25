      SUBROUTINE colloid_compute_interaction_implicit_velocity_pair_solver(this,&
           comm,MPI_PREC,dt,drag,torque,FB,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_interaction_implicit_
        !               velocity_pair_solver
        !----------------------------------------------------
        !
        ! Purpose     : Update velocties using implicit
        !               splitting scheme for 
        !               pair-wise colloids, using lapack.
        !               Update positions using explicit 
        !               scheme.
        !               
        !
        !               The interactions may include:
        !               1)lubrication force correction,
        !               2)repulsive force to prevent
        !               overlap(contact force).
        !               One can also use 2) as physical
        !               DLVO repulsive force.
        !
        ! Routines    :
        !
        ! References  : Litvinov et al. J. Comput. Phys. 2010.
        !
        ! Remarks     : Consider only serial run.
        !
        ! Revisions   : V0.1 09.05.2012, original version.
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
        REAL(MK), INTENT(IN)                    :: dt
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: drag
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: torque
        REAL(MK), DIMENSION(:,:), INTENT(OUT)   :: FB
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        ! dim : dimension
        ! num : number of total colloids
        ! ent : number of unknown entries in left hand matrix.
        ! ML  : left hand matrix
        ! MR  : right hand matrix
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim, num, ent, dim2        
        INTEGER                                 :: num_sweep, s
        REAL(MK)                                :: dt_s
        REAL(MK)                                :: hn_l, hm_l
        REAL(MK)                                :: hn_r, hm_r
        REAL(MK)                                :: F0, F1
        REAL(MK)                                :: F0_repul

        INTEGER                                 :: i, j, k, m        
        
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: ML, MR
        INTEGER, DIMENSION(:), ALLOCATABLE      :: IPIV
        REAL(MK)                                :: ai, aj, aa, Aij, Bij 
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: F_lub
        INTEGER                                 :: num_wall_solid

        
        REAL(MK),DIMENSION(3)                   :: x_image, v_image        
        REAL(MK), DIMENSION(3)                  :: rij, eij, vij
        REAL(MK)                                :: ve
        REAL(MK)                                :: r, h
        REAL(MK),DIMENSION(3)                   :: Fi_repul
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: F_repul
        REAL(MK),DIMENSION(3,6)                 :: FB_lub
        REAL(MK),DIMENSION(3,6)                 :: FB_repul
      
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------

        stat_info     = 0
        stat_info_sub = 0
        
        dim   = this%num_dim
        num   = this%num_colloid
        ent   = dim * 2
        dim2  = dim * 2

        num_sweep = this%implicit_pair_num_sweep
        dt_s      = dt / num_sweep
    
        hn_l  = this%cc_lub_cut_off
        hm_l  = this%cc_lub_cut_on
        hn_r  = this%cc_repul_cut_off
        hm_r  = this%cc_repul_cut_on
        F0_repul = this%cc_repul_F0
        
        !----------------------------------------------------
        ! For now, we assume that each radius is the same.
        !----------------------------------------------------
        
        ai = this%radius(1,1)
        aj = ai
        aa = ai+aj
        
        IF ( dim == 2 ) THEN
           
           F0  = 3.0_MK*mcf_pi*SQRT(2.0_MK)/4.0_MK
           F1  = 231.0_MK*mcf_pi*SQRT(2.0_MK) / 80.0_MK
           
        ELSE
           
           PRINT *, __FILE__, __LINE__, "dimension not available!"
           stat_info = -1
           GOTO 9999
           
        END IF
         
        ALLOCATE(F_lub(dim,num))
       
        ALLOCATE(ML(1:ent, 1:ent))
        ALLOCATE(MR(1:ent, 1))
        ALLOCATE(IPIV(1:ent))
        ML(:,:) = 0.0_MK
        MR(:,:) = 0.0_MK
        IPIV(:) = 0
        
        num_wall_solid = &
             boundary_get_num_wall_solid(this%boundary,stat_info_sub)
        
        ALLOCATE(F_repul(dim,num))
        
        
        !----------------------------------------------------
        ! Do the pairwise sweep num_sweep times to update
        ! velocities using only lubrcation correction forces.
        ! Time step at each sweep is dt_s=dt/num_sweep.
        ! First half dt_s is explicit and second
        ! half dt_s is implicit.
        !
        ! Also inside the loop record repulsive forces
        ! between colloid-colloid and colloid-wall, which
        ! are used later for updating velocity explicitly.
        !----------------------------------------------------
        
        DO s = 1, num_sweep
           
           !-------------------------------------------------
           ! Reset lubrication correction forces at beginning
           ! of each sweeping.
           !-------------------------------------------------
           
           F_lub(1:dim,1:num)   = 0.0_MK
                   
           !-------------------------------------------------
           ! Explicit: first half dt_s.
           !
           ! Loop over each colloid and its interaction with 
           ! other colloids and 
           ! get lubrication correction forces.
           !-------------------------------------------------
           
           DO i = 1, num - 1
              
              DO j = i + 1, num
                 
                 !-------------------------------------------
                 ! Calculate the gap and unit vector joining
                 ! two colloids, images of colloids have to
                 ! be considered.
                 !-------------------------------------------
                 
                 CALL colloid_nearest_image(this,&
                      this%x(1:dim,i),j, &
                      x_image(1:dim),rij(1:dim), &
                      v_image(1:dim),stat_info_sub)
                 
                 r  = SQRT(DOT_PRODUCT(rij(1:dim), rij(1:dim)))
                 h  = r - aa
                 eij(1:dim) = rij(1:dim) / r
                 vij(1:dim) = this%v(1:dim,i,1) - v_image(1:dim)
                 ve = DOT_PRODUCT(vij(1:dim), eij(1:dim))
                 
                 !-------------------------------------------
                 ! If gap is smaller than hn_l, i.e., cut_off
                 ! of lubrication correction,
                 ! it needs lubrication correction.
                 !-------------------------------------------
                 
                 IF ( h < hn_l ) THEN
                    
                    !----------------------------------------
                    ! If gap is smaller than minimal allowed
                    ! gap, set it to the pre-set minimum.
                    !----------------------------------------
                    
                    IF ( h < hm_l ) THEN
                       
                       h = hm_l
                       
                    END IF
                    
                    Aij = -0.5_MK * this%eta * ve * &
                         ( (aa/h)**1.5_MK  * (F0 + h*F1/aa) -&
                         (aa/hn_l)**1.5_MK * (F0 + hn_l*F1/aa) )
                    
                    F_lub(1:dim,i) = F_lub(1:dim,i) + Aij * eij(1:dim)
                    F_lub(1:dim,j) = F_lub(1:dim,i) - Aij * eij(1:dim)
                    
                 END IF ! h < hn_l

              END DO ! j = i + 1, num
              
           END DO ! i = 1, num - 1
              
           !-------------------------------------------------
           ! After knowing lubrication correction forces
           ! for all colloids, integrate velocities using 
           ! explicit scheme by dt_s/2.
           !-------------------------------------------------
           
           DO i = 1, num
              
              this%v(1:dim,i,1) = this%v(1:dim,i,1) + &
                   0.5_MK * F_lub(1:dim,i) * dt_s / this%m(i)
              
           END DO ! i = 1, num
           
           !-------------------------------------------------
           ! Implicit: second half dt_s.
           !          
           ! Loop over each colloid and its interaction with 
           ! other colloids.
           ! Build up left hand matrix and right hand vector
           ! for a system of linear equations of two collodis,
           ! and then update velocity using 
           ! implicit scheme by dt_s/2.
           !-------------------------------------------------
           
           DO i = 1, num - 1
              
              !----------------------------------------------
              ! Construct the left hand side in the two loops.
              !----------------------------------------------
              
              !----------------------------------------------
              ! Self interaction.
              !----------------------------------------------
              
              DO k = 1, ent
                 
                 ML(k,k) = 1.0_MK
                 
              END DO
              
              !----------------------------------------------
              ! Interaction between two different colloids.
              !----------------------------------------------
            
              DO j = i + 1, num
                 
                 !-------------------------------------------
                 ! Construct right hand side from 
                 ! available current velocities.
                 !-------------------------------------------
                 
                 MR(1:dim,1)       = this%v(1:dim,i,1)
                 MR(dim+1:2*dim,1) = this%v(1:dim,j,1)
                 
                 !-------------------------------------------
                 ! Calculate the gap and unit vector joining
                 ! two colloids, images of colloids have to
                 ! be considered according to different
                 ! boundary conditions.
                 !-------------------------------------------
                 
                 CALL colloid_nearest_image(this,&
                      this%x(1:dim,i),j, &
                      x_image(1:dim),rij(1:dim), &
                      v_image(1:dim),stat_info_sub)
                 
                 r  = SQRT(DOT_PRODUCT(rij(1:dim), rij(1:dim)))
                 h  = r - aa
                 eij(1:dim) = rij(1:dim) / r
              
                 !-------------------------------------------
                 ! If gap is smaller than hn_l, i.e., cut_off
                 ! of lubrication correction,
                 ! it needs lubrication correction.
                 !-------------------------------------------
              
                 IF ( h < hn_l ) THEN
                    
                    !----------------------------------------
                    ! If gap is smaller than minimal allowed
                    ! gap, set it to the pre-set minimum.
                    !----------------------------------------
                 
                    IF ( h < hm_l ) THEN
                       
                       h = hm_l
                       
                    END IF
                    
                    !----------------------------------------
                    ! dt_s/2 is coefficient of Aij.
                    !----------------------------------------
                    
                    Aij = -0.5_MK * this%eta * &
                         ( (aa/h)**1.5_MK  * (F0 + h*F1/aa) -&
                         (aa/hn_l)**1.5_MK * (F0 + hn_l*F1/aa) ) * &
                         0.5_MK * dt_s
                    
                    !----------------------------------------
                    ! Divided by mass of colloid to have
                    ! correct unit of velocity.
                    !----------------------------------------

                    DO k = 1, dim
                       
                       DO m = 1, dim
                          
                          ML(k,m) = &
                               ML(k,m) - &
                               Aij * eij(k) * eij(m) / this%m(i)
                          
                          ML(k,dim+m) = &
                               ML(k,dim+m) + &
                               Aij * eij(k) * eij(m) / this%m(i)
                          
                       END DO ! m
                       
                    END DO ! k
                    
                    DO k = 1, dim
                       
                       DO m = 1, dim
                          
                          ML(dim+k,dim+m) = &
                               ML(dim+k,dim+m) - &
                               Aij * eij(k) * eij(m) / this%m(j)
                          
                          ML(dim+k,m) = &
                               ML(dim+k,m) + &
                               Aij * eij(k) * eij(m) / this%m(j)
                          
                       END DO ! m
                       
                    END DO ! k
                    
                    !----------------------------------------
                    ! Solving system of linear equations
                    ! L*x=R,
                    ! where x is the solution for v^(n+1).
                    !----------------------------------------
                    
                    !PRINT *, "ML, MR: ", ML(:,:), MR(:,1)
                    CALL tool_solve_linear_equations(this%tool, &
                         ent,1, ML, ent, IPIV, MR, ent,stat_info_sub)
                    !PRINT *, "MR: ", MR(:,1)
                    
                    IF ( stat_info_sub /=0 ) THEN
                       
                       PRINT *, "ML, MR: ", ML(:,:), MR(:,1)
                       PRINT *, __FILE__, __LINE__, &
                            "solving linear equations failed!"
                       stat_info = -1
                       GOTO 9999
                       
                    END IF
                 
                    !----------------------------------------
                    ! update the velocity of the interacting
                    ! pair immediately.
                    !----------------------------------------
                    
                    this%v(1:dim,i,1) = MR(1:dim,1)
                    this%v(1:dim,j,1) = MR(dim+1:2*dim,1)
                    
                 END IF ! h < hn_l
                 
              END DO ! j = 1, num
              
           END DO ! i = 1, num -1
           
        END DO ! s = 1, num_sweep
        
        !STOP
        !----------------------------------------------------
        ! Loop over each colloid and its interaction with 
        ! other colloids and walls,
        ! calculate the repulsive force between 
        ! colloid-colloid and colloid-wall.
        !
        ! Update velocity explicitly from 
        ! 1: SPH solvent forces
        ! 2: pairwise repulsive force between colloid-colloid
        ! 3: repulsive force between colloid-wall
        ! by two steps, each one dt/2.
        !----------------------------------------------------
        
        
        !----------------------------------------------------
        ! Repulsive forces between colloid-colloid
        !----------------------------------------------------

        F_repul(1:dim,1:num) = 0.0_MK

        DO i = 1, num - 1
           
           !-------------------------------------------------
           ! Interaction between two different colloids.
           !-------------------------------------------------
           
           DO j = i + 1, num
              
              !----------------------------------------------
              ! Calculate the gap and unit vector joining
              ! two colloids, images of colloids have to
              ! be considered according to different
              ! boundary conditions.
              !----------------------------------------------
                 
              CALL colloid_nearest_image(this,&
                   this%x(1:dim,i),j, &
                   x_image(1:dim),rij(1:dim), &
                   v_image(1:dim),stat_info_sub)
              
              r  = SQRT(DOT_PRODUCT(rij(1:dim), rij(1:dim)))
              h  = r - aa
              eij(1:dim) = rij(1:dim) / r
        
              !----------------------------------------------
              ! If gap is smaller than hn_r, i.e., cut_off 
              ! of repulsive force, it needs repulsive force.
              !----------------------------------------------
                 
              IF ( h < 5.0_MK * hn_r ) THEN
                 
                 Bij = 0.0_MK
                 
                 SELECT CASE ( this%cc_repul_type )
                    
                 CASE ( mcf_cc_repul_type_Hookean )
                    !----------------------------------------
                    ! For linear spring force, it has clear 
                    ! zero at hn_r.
                    !----------------------------------------
                    
                    IF ( h < hn_r ) THEN
                       
                       !-------------------------------------
                       ! If gap smaller than minimal allowed 
                       ! gap, set it to the pre-set minimum.
                       !-------------------------------------
                       
                       IF ( h < hm_r ) THEN
                          
                          h = hm_r
                          
                       END IF
                       
                       Bij = F0_repul - F0_repul*h/hn_r
                       
                    END IF ! h < hn_r
                    
                 CASE ( mcf_cc_repul_type_DLVO )
                    
                    !----------------------------------------
                    ! For DLVO force, it does not have clear 
                    ! zero at hn. But at 5*hn, its value 
                    ! is smaller than F0/100.
                    !----------------------------------------
                       
                    !----------------------------------------
                    ! If gap smaller than minimal allowed 
                    ! gap, set it to the pre-set minimum.
                    !----------------------------------------
                       
                    IF ( h < hm_r ) THEN
                       
                       h = hm_r
                       
                    END IF
                    
                    Bij = F0_repul / hn_r * &
                         EXP(-h/hn_r) /(1.0_MK-EXP(-h/hn_r))
                    
                 END SELECT ! cc_repul_type
                 
                 !-------------------------------------------
                 ! Save pair-wise repulsive forces
                 !-------------------------------------------
                 
                 F_repul(1:dim,i) = F_repul(1:dim,i) + &
                      Bij * eij(1:dim)
                 F_repul(1:dim,j) = F_repul(1:dim,j) - &
                      Bij * eij(1:dim)
                 
              END IF ! h < hn_r
              
           END DO ! j = 1, num
           
        END DO ! i = 1, num -1
        
              
        !----------------------------------------------------
        ! Repulsive forces between colloid-wall
        !----------------------------------------------------
        
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
                   F_repul(1:dim, i ) + Fi_repul(1:dim)
              
           END DO ! i = 1, num
           
        END IF ! cw_repul_type > no
        
        !----------------------------------------------------
        ! Update velocity using SPH solvent force and 
        ! repulsive force by dt/2
        !----------------------------------------------------
        
        DO i = 1, num
           
           this%v(1:dim,i,1) = &
                this%v(1:dim,i,1) + &
                0.5_MK * ( F_repul(1:dim,i)+drag(1:dim,i) ) * dt /&
                this%m(i)           
           
        END DO ! i = 1, num

        !----------------------------------------------------
        ! Update position by dt
        !----------------------------------------------------
        
        DO i = 1, num
           
           this%x(1:dim,i) = this%x(1:dim,i) + &
                this%v(1:dim,i,1) * dt
           
        END DO
        
        !----------------------------------------------------
        ! Calculate the pairwise repulsive force again 
        ! using new positions.
        !----------------------------------------------------
        
        F_repul(1:dim,1:num)   = 0.0_MK
        
        DO i = 1, num - 1
           
           DO j = i + 1, num
              
              !----------------------------------------------
              ! Calculate the gap and unit vector joining
              ! two colloids, images of colloids have to
              ! be considered according to different
              ! boundary conditions.
              !----------------------------------------------
                 
              CALL colloid_nearest_image(this,&
                   this%x(1:dim,i),j, &
                   x_image(1:dim),rij(1:dim), &
                   v_image(1:dim),stat_info_sub)
              
              r  = SQRT(DOT_PRODUCT(rij(1:dim), rij(1:dim)))
              h  = r - aa
              eij(1:dim) = rij(1:dim) / r
              
              !----------------------------------------------
              ! If gap is smaller than hn_r, i.e., cut_off 
              ! of repulsive force.
              ! it needs repulsive force.
              !----------------------------------------------
              
              IF ( h < 5.0_MK * hn_r ) THEN
                 
                 Bij = 0.0_MK
                 
                 SELECT CASE ( this%cc_repul_type )
                    
                 CASE ( mcf_cc_repul_type_Hookean )
                    !----------------------------------------
                    ! For linear spring force, it has clear 
                    ! zero at hn_r.
                    !----------------------------------------
                    
                    IF ( h < hn_r ) THEN
                       
                       !----------------------------------
                       ! If gap smaller than minimal allowed 
                       ! gap, set it to the pre-set minimum.
                       !----------------------------------
                       
                       IF ( h < hm_r ) THEN
                          
                          h = hm_r
                          
                       END IF
                       
                       Bij = F0_repul - F0_repul*h/hn_r
                       
                    END IF ! h < hn_r
                    
                 CASE ( mcf_cc_repul_type_DLVO )
                    
                    !-------------------------------------
                    ! For DLVO force, it does not have clear 
                    ! zero at hn. But at 5*hn, its value 
                    ! is smaller than F0/100.
                    !-------------------------------------
                    
                    !-------------------------------------
                    ! If gap smaller than minimal allowed 
                    ! gap, set it to the pre-set minimum.
                    !-------------------------------------
                    
                    IF ( h < hm_r ) THEN
                       
                       h = hm_r
                       
                    END IF
                    
                    Bij = F0_repul / hn_r * &
                         EXP(-h/hn_r) /(1.0_MK-EXP(-h/hn_r))
                    
                 END SELECT ! cc_repul_type
                 
                 
                 !----------------------------------------
                 ! Save pair-wise repulsive forces
                 !----------------------------------------
                 
                 F_repul(1:dim,i) = F_repul(1:dim,i) + &
                      Bij * eij(1:dim)
                 F_repul(1:dim,j) = F_repul(1:dim,j) - &
                      Bij * eij(1:dim)
                 
              END IF ! h < hn_r
              
           END DO ! j = i+1, num
           
        END DO ! i = 1, num-1
        
        !----------------------------------------------------
        ! Repulsive force between colloid-wall.
        ! Save the repulsive force on the wall also.
        !----------------------------------------------------
        
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
                   F_repul(1:dim, i ) + Fi_repul(1:dim)
              FB(1:dim,1:dim2) = &
                   FB(1:dim,1:dim2) + FB_repul(1:dim,1:dim2)
              
           END DO ! i = 1, num
           
        END IF ! cw_repul_type > no
        
        !----------------------------------------------------
        ! Update velocity using SPH solvent force and 
        ! repulsive force at second dt/2
        !----------------------------------------------------
        
        DO i = 1, num
           
           this%v(1:dim,i,1) = &
                this%v(1:dim,i,1) + &
                0.5_MK * ( F_repul(1:dim,i)+drag(1:dim,i) ) * dt /&
                this%m(i)           
           
        END DO ! i = 1, num

        !----------------------------------------------------
        ! Update total torque on each colloid
        !----------------------------------------------------
        
        this%torque(1:3,1:num) = torque(1:3,1:num)
        
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE colloid_compute_interaction_implicit_velocity_pair_solver
      
      
      
