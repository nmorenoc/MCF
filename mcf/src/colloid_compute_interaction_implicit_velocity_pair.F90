      SUBROUTINE colloid_compute_interaction_implicit_velocity_pair(this,&
           comm,MPI_PREC,dt,drag,torque,FB,stat_info)
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
        ! Revisions   : V0.2 18.05.2012, shift 2nd repulsive force
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
        INTEGER                                 :: dim, num, dim2
        
        REAL(MK)                                :: hn_l, hm_l
        REAL(MK)                                :: F0, F1
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
        
        INTEGER                                 :: num_sweep, sweep
        REAL(MK)                                :: dt_sweep
        
        INTEGER                                 :: num_sub_step, sub_step
        REAL(MK)                                :: dt_sub
        
        REAL(MK)                                :: ai, aj, aa
        REAL(MK)                                :: Aij
        REAL(MK), DIMENSION(3)                  :: vi_old, vj_old
        REAL(MK), DIMENSION(3)                  :: vi_new, vj_new
        REAL(MK)                                :: Bij
        
        REAL(MK),DIMENSION(3)                   :: x_image, v_image 
        REAL(MK), DIMENSION(3)                  :: rij, eij, vij
        REAL(MK)                                :: r, h, ve
          
        REAL(MK)                                :: fa, dt_f
        INTEGER                                 :: i, j, k, m        
        
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------

        stat_info     = 0
        stat_info_sub = 0
        
        FB(:,:) = 0.0_MK

        !----------------------------------------------------
        ! Back up previous velocities.
        !----------------------------------------------------

        ALLOCATE(v_backup(dim,num))
        v_backup(1:dim,1:num) = this%v(1:dim,1:num,1)
    
        !----------------------------------------------------
        ! set up physics parameters.
        !----------------------------------------------------

        dim   = this%num_dim
        num   = this%num_colloid
        dim2  = dim * 2
        
        hn_l  = this%cc_lub_cut_off
        hm_l  = this%cc_lub_cut_on
        
        IF ( dim == 2 ) THEN
           
           F0  = 3.0_MK*mcf_pi*SQRT(2.0_MK)/4.0_MK
           F1  = 231.0_MK*mcf_pi*SQRT(2.0_MK) / 80.0_MK
           
        ELSE
           
           PRINT *, __FILE__, __LINE__, "dimension not available!"
           stat_info = -1
           GOTO 9999
           
        END IF
        
        hn_r  = this%cc_repul_cut_off
        hm_r  = this%cc_repul_cut_on
        F0_repul = this%cc_repul_F0

        ALLOCATE(F_repul(dim,num))
        F_repul(:,:) = 0.0_MK
  
        num_wall_solid = &
             boundary_get_num_wall_solid(this%boundary,stat_info_sub)
     
        !----------------------------------------------------
        ! set up numerical parameters.
        !----------------------------------------------------
        
        num_sweep = this%implicit_pair_num_sweep
        dt_sweep  = dt / num_sweep
        
        num_sub_step = this%explicit_sub_time_step
        dt_sub       = dt / num_sub_step
        
        !----------------------------------------------------
        ! For now, we assume that each radius is the same.
        !----------------------------------------------------
        
        ai = this%radius(1,1)
        aj = ai
        aa = ai+aj
        
        !----------------------------------------------------
        ! First update velocity from contribution SPH forces.
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
        ! Do the pairwise sweep num_sweep times to update
        ! velocities using only lubrcation correction forces.
        ! Time step at each sweep is dt_sweep=dt/num_sweep.
        !
        ! It is implicit.
        !
        !----------------------------------------------------
        
        IF ( this%cc_lub_type > mcf_cc_lub_type_no ) THEN
           
           DO sweep = 1, num_sweep
              
              !----------------------------------------------
              ! Implicit: dt_sweep.
              !          
              ! Loop over each colloid and its interaction with 
              ! other colloids.
              ! Build up left hand matrix and right hand vector
              ! for a system of linear equations of two collodis,
              ! and then update velocity using implicit scheme
              ! by dt_sweep.
              ! We solve this small matrix explicitly/analytically
              ! using maxima.
              !----------------------------------------------
              
              DO i = 1, num - 1
                 
                 !-------------------------------------------
                 ! Interaction between two different colloids.
                 !-------------------------------------------
            
                 DO j = i + 1, num
                    
                    vi_old(1:dim) = this%v(1:dim,i,1)
                    vj_old(1:dim) = this%v(1:dim,j,1)
                    
                    !----------------------------------------
                    ! Calculate the gap and unit vector joining
                    ! two colloids, images of colloids have to
                    ! be considered according to different
                    ! boundary conditions.
                    !----------------------------------------
                    
                    CALL colloid_nearest_image(this,&
                         this%x(1:dim,i),j, &
                         x_image(1:dim),rij(1:dim), &
                         v_image(1:dim),stat_info_sub)
                    
                    r  = SQRT(DOT_PRODUCT(rij(1:dim), rij(1:dim)))
                    h  = r - aa
                    eij(1:dim) = rij(1:dim) / r
                    
                    !----------------------------------------
                    ! If gap is smaller than hn_l, i.e., 
                    ! cut_off of lubrication correction,
                    ! it needs lubrication correction.
                    !----------------------------------------
                    
                    IF ( h < hn_l ) THEN
                       
                       !-------------------------------------
                       ! If gap is smaller than minimal allowed
                       ! gap, set it to the pre-set minimum.
                       !-------------------------------------
                 
                       IF ( h < hm_l ) THEN
                          
                          h = hm_l
                          
                       END IF
                       
                       !-------------------------------------
                       ! dt_sweep/2/m is coefficient of Aij.
                       ! Assuming now mi=mj
                       !-------------------------------------
                       
                       Aij = -0.5_MK * this%eta * &
                            ( (aa/h)**1.5_MK  * (F0 + h*F1/aa) - &
                            (aa/hn_l)**1.5_MK * (F0 + hn_l*F1/aa) ) * &
                            dt_sweep / this%m(i)
                       
                       !-------------------------------------
                       ! Solve 2*2 linear system analytically.
                       !-------------------------------------
                       
#include "velocity_pair_backward_Euler_2d.inc"
                       
                       !-------------------------------------
                       ! update velocity of the interacting
                       ! pair immediately.
                       !-------------------------------------
                       
                       this%v(1:dim,i,1) = vi_new(1:dim)
                       this%v(1:dim,j,1) = vj_new(1:dim)
                       
                    END IF ! h < hn_l
                    
                 END DO ! j = i+1, num
                 
              END DO ! i = 1, num-1
              
           END DO ! s = 1, num_sweep
           
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
        ! Repulsive forces between colloid-colloid and
        ! colloid-wall.
        !----------------------------------------------------
        
        IF ( this%cc_repul_type > mcf_cc_repul_type_no .OR. &
             this%cw_repul_type > mcf_cw_repul_type_no ) THEN
           
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
                   sub_step == 1) THEN
                 
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
                   sub_step == 1) THEN
                 
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
                       
                    END DO ! j = i+1, num
                    
                 END DO ! i = 1, num-1
                 
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
           
        END IF ! cc_repul_type > no AND cw_repul_type > no
        
        !----------------------------------------------------
        ! Work out the effective drag by calculating the
        ! difference between new and old velocity and
        ! assuming an Euler scheme was used as integrator.
        !----------------------------------------------------
        
        DO i = 1, num
           
           this%drag(1:dim,i) = &
                (this%v(1:dim,i,1) - v_backup(1:dim,i))*this%m(i) / dt
           
        END DO
        
        !----------------------------------------------------
        ! Update total torque on each colloid
        !----------------------------------------------------
        
        this%torque(1:3,1:num) = torque(1:3,1:num)
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE colloid_compute_interaction_implicit_velocity_pair
      
      
      
