      SUBROUTINE colloid_compute_interaction_implicit_velocity_all(this,&
           comm,MPI_PREC,dt,drag,torque,FB,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_interaction_implicit_
        !               velocity_all
        !----------------------------------------------------
        !
        ! Purpose     : Calculate updated velocties using 
        !               old velocities by implicit scheme for 
        !               all colloids at once.
        !               
        !               Afterwards, compute colloid-colloid
        !               interactions using updated velocities.
        !               In the end, compute colloid-wall
        !               interactions.
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
        ! References  :
        !
        ! Remarks     : Consider only serial run.
        !
        ! Revisions   : V0.1 07.05.2012, original version.
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
        ! ent : number of unknown entries in the left hand matrix.
        ! ML  : left hand matrix
        ! MR  : right hand matrix
        !----------------------------------------------------
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: i, j, k, m
        INTEGER                                 :: dim, num, ent
        INTEGER                                 :: dim2
        REAL(MK),DIMENSION(3)                   :: x_image, v_image        
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: ML, MR
        REAL(MK), DIMENSION(3)                  :: rij, eij, vij
        REAL(MK)                                :: r, h
        REAL(MK)                                :: hn_l, hm_l
        REAL(MK)                                :: hn_r, hm_r
        LOGICAL                                 :: interact_l, interact_r
        REAL(MK)                                :: F0, F1, F_rep
        REAL(MK)                                :: ai, aj, aa, Aij, Bij
        INTEGER, DIMENSION(:), ALLOCATABLE      :: IPIV
        
        REAL(MK),DIMENSION(3)                   :: F_i
        REAL(MK),DIMENSION(3,6)                 :: FB_lub
        REAL(MK),DIMENSION(3,6)                 :: FB_repul
        INTEGER                                 :: num_wall_solid
      
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------

        stat_info = 0
        stat_info_sub = 0
                
        dim = this%num_dim
        num = this%num_colloid
        ent = dim * num
        dim2 = dim * 2
        hn_l  = this%cc_lub_cut_off
        hm_l  = this%cc_lub_cut_on
        hn_r  = this%cc_repul_cut_off
        hm_r  = this%cc_repul_cut_on
        F_rep = this%cc_repul_F0
        
        ALLOCATE(ML(1:ent, 1:ent))
        ALLOCATE(MR(1:ent, 1))
        ALLOCATE(IPIV(1:ent))
        ML(:,:) = 0.0_MK
        MR(:,:) = 0.0_MK
        IPIV(:) = 0
        
        num_wall_solid = &
             boundary_get_num_wall_solid(this%boundary,stat_info_sub)
        
        FB(:,:)       = 0.0_MK
        
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
        
        !----------------------------------------------------
        ! Loop over each colloid and its interaction with 
        ! other colloids and walls.
        ! Build up the left hand matrix and right hand vector
        ! for the system of linear equations
        !----------------------------------------------------
        
        DO i = 1, num
           
           !-------------------------------------------------
           ! Construct the right hand side
           ! 1st contribution is the old velocity
           ! 2nd contribution is from the SPH solvent force.
           !-------------------------------------------------
           
           MR((i-1)*dim+1:i*dim,1) = &
                this%v(1:dim,i,1) + drag(1:dim,i) * dt / this%m(i)
           
           !-------------------------------------------------
           ! Construct the left hand side in the two loops.
           !-------------------------------------------------

           DO j = 1, num
              
              !----------------------------------------------
              ! Two colloids are the same one.
              !----------------------------------------------
              
              IF ( j == i ) THEN
                 
                 DO k = 1, dim
                    
                    ML((i-1)*dim+k,(i-1)*dim+k) = 1.0_MK
                    
                 END DO
                 
                 !----------------------------------------------
                 ! Two different colloids
                 !----------------------------------------------
                 
              ELSE

                 
                 !----------------------------------------------
                 ! Calculate the gap and unit vector joining
                 ! two colloids, images of colloids have to
                 ! be considered.
                 !----------------------------------------------

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
                 
                    Aij = -0.5_MK * this%eta * &
                         ( (aa/h)**1.5_MK  * (F0 + h*F1/aa) -&
                         (aa/hn_l)**1.5_MK * (F0 + hn_l*F1/aa) ) * &
                         dt / this%m(i)
                    
                    
                    DO k = 1, dim
                       
                       DO m = 1, dim
                          
                          ML((i-1)*dim+k,(i-1)*dim+m) = &
                               ML((i-1)*dim+k,(i-1)*dim+m) - &
                               Aij * eij(k) * eij(m)
                          
                          ML((i-1)*dim+k,(j-1)*dim+m) = &
                               ML((i-1)*dim+k,(j-1)*dim+m) + &
                               Aij * eij(k) * eij(m)
                          
                       END DO ! m
                       
                    END DO ! k
                    
                 END IF ! h < hn_l
                 
                 !-------------------------------------------
                 ! If gap is smaller than hn_r, i.e., cut_off 
                 ! of repulsive force.
                 ! it needs repulsive force.
                 !-------------------------------------------
                 
                 IF ( h < 5.0_MK * hn_r ) THEN
                    
                    SELECT CASE ( this%cc_repul_type )
                       
                    CASE ( mcf_cc_repul_type_Hookean )
                       !-------------------------------------
                       ! For linear spring force, it has clear 
                       ! zero at hn_r.
                       !-------------------------------------
                       
                       IF ( h < hn_r ) THEN
                          
                          !----------------------------------
                          ! If gap smaller than minimal allowed 
                          ! gap, set it to the pre-set minimum.
                          !----------------------------------
                          
                          IF ( h < hm_r ) THEN
                             
                             h = hm_r
                             
                          END IF
                          
                          Bij = F_rep - F_rep*h/hn_r
                          
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
                       
                       Bij = F_rep / hn_r * &
                            EXP(-h/hn_r) /(1.0_MK-EXP(-h/hn_r))
                       
                    END SELECT
                    
                    !----------------------------------------
                    ! Construct the right hand side
                    ! 3rd contribution is from the repulsion
                    ! between colloid-colloid
                    !----------------------------------------
        
                    MR((i-1)*dim+1:i*dim,1) = &
                         MR((i-1)*dim+1:i*dim,1) + &
                         Bij*dt*eij(1:dim)/this%m(i)
                    
                 END IF ! h < hn_r
              
              END IF ! j == i and j /= i
              
           END DO ! j = 1, num
           
           !-------------------------------------------------
           ! Compute interactions between colloid-wall and
           ! put it in right hand side.
           !-------------------------------------------------
           
           IF ( num_wall_solid > 0 .AND. &
                ( this%cw_lub_type > mcf_cw_lub_type_no .OR. &
                this%cw_repul_type > mcf_cw_repul_type_no) ) THEN
              
              !----------------------------------------------
              ! If lubrication correction is required to 
              ! amend forces between colloid and wall, 
              ! forces are introduced for colloid-wall gap
              ! smaller than cw_lub_cut_off.
              !
              ! 2D, 3D lurication theory formulations
              ! are different.
              !----------------------------------------------
              
              FB_lub(:,:)   = 0.0_MK
              FB_repul(:,:) = 0.0_MK
              
              IF ( this%cw_lub_type > mcf_cw_lub_type_no ) THEN
                 
                 CALL colloid_compute_lubrication_cw(this,&
                      this%x(1:dim,i),this%v(1:dim,i,1),i,&
                      F_i(1:dim),FB_lub(1:dim,1:dim2),stat_info_sub)
                 
                 IF ( stat_info_sub /= 0 ) THEN
                    PRINT *, __FILE__, __LINE__, &
                         "calling lubrication cw failed!"
                    stat_info = -1
                    GOTO 9999
                 END IF

                 !-------------------------------------------
                 ! Construct the right hand side
                 ! 4th contribution is from the lubrication
                 ! correction between colloid-wall.
                 !-------------------------------------------
                 
                 MR((i-1)*dim+1:i*dim,1) = &
                      MR((i-1)*dim+1:i*dim,1) + &
                      dt*F_i(1:dim)/this%m(i)
                 
              END IF ! cw_lub_type
              
              !----------------------------------------------
              ! Add up repulsive force to prevent overlaps
              ! between wall and colloid.
              !----------------------------------------------
              
              IF ( this%cw_repul_type > mcf_cw_repul_type_no ) THEN
                 
                 CALL colloid_compute_repulsion_cw(this,&
                      this%x(1:dim,i),i,F_i(1:dim),&
                      FB_repul(1:dim,1:dim2),stat_info_sub)
                 
                 IF ( stat_info_sub /= 0 ) THEN
                    PRINT *, __FILE__, __LINE__, &
                         "calling repulsion cw failed !"
                    stat_info = -1
                    GOTO 9999
                 END IF
        
                 !-------------------------------------------
                 ! Construct the right hand side
                 ! 5th contribution is from the repulsion
                 ! between colloid-wall.
                 !-------------------------------------------
                 
                 MR((i-1)*dim+1:i*dim,1) = &
                      MR((i-1)*dim+1:i*dim,1) + &
                      dt*F_i(1:dim)/this%m(i)
                 
              END IF ! cc_repul_type              
              
              FB(1:dim,1:dim2) =  FB(1:dim,1:dim2) + &
                   FB_lub(1:dim,1:dim2) + FB_repul(1:dim,1:dim2)
              
           END IF ! num_wall_solid > 0
           
        END DO ! i = 1, num
        
        !----------------------------------------------------
        ! Solving system of linear equations
        ! L*x=R,
        ! where x is the solution for v^(n+1).
        !----------------------------------------------------
        
        PRINT *, "ML, MR: ", ML(:,:), MR(:,1)
        CALL tool_solve_linear_equations(this%tool, &
             ent,1, ML, ent, IPIV, MR, ent,stat_info_sub)
        PRINT *, "MR: ", MR(:,1)
        !STOP
        IF ( stat_info_sub /=0 ) THEN
           
           PRINT *, "ML, MR: ", ML(:,:), MR(:,1)
           PRINT *, __FILE__, __LINE__, &
                "solving linear equations failed!"
           stat_info = -1
           GOTO 9999
           
        END IF

        !----------------------------------------------------
        ! Using updated velocity and old velocity to calcuate
        ! the effective total force at this time step.
        !----------------------------------------------------
        
        DO i = 1, num
           
           DO j = 1, dim
              
              k = (i-1)*dim+j
              this%drag(j,i) = &
                   this%m(i) * ( MR(k,1) - this%v(j,i,1)) / dt
              
           END DO
           
        END DO
        
        !----------------------------------------------------
        ! Update total torque on each colloid
        !----------------------------------------------------

        this%torque(1:3,1:num) = torque(1:3,1:num)
        
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE colloid_compute_interaction_implicit_velocity_all
      
      
      
