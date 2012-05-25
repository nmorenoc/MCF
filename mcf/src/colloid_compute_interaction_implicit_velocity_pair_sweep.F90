      SUBROUTINE colloid_compute_interaction_implicit_velocity_pair_sweep(&
           this, dt, num_sweep, stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_interaction_implicit_
        !               velocity_pair_sweep
        !----------------------------------------------------
        !
        ! Purpose     :   Update velocties from 
        !                 lubrication-correction force using 
        !                 implicit splitting scheme for 
        !                 pair-wise colloids with number of
        !                 sweeps.
        !         
        !
        ! Routines    :
        !
        ! References  : Shardlow T. SIAM J. Sci. Comput. 2003.
        !               Litvinov S. et al. J. Comput. Phys. 2010.
        !
        ! Remarks     :
        !
        ! Revisions   : V0.1 25.05.2012, original version.
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
        REAL(MK), INTENT(IN)                    :: dt
        INTEGER, INTENT(IN)                     :: num_sweep
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local physics variables:
        !
        ! dim   : dimension
        ! num   : number of total colloids
        !
        ! hn_l,hm_l : lubrication correction cut off and on.
        ! F0,F1     : lubrication correction parameters.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim, num
        REAL(MK)                                :: hn_l, hm_l
        REAL(MK)                                :: F0, F1
        
        !----------------------------------------------------
        ! Local numerical variables:
        !
        ! num_sweep   : number of sweeps for pair-wise implicit
        !               splitting scheme.
        !----------------------------------------------------
        INTEGER                                 :: sweep
        REAL(MK)                                :: dt_sweep
        
        !----------------------------------------------------
        ! ai,aj       : radius of colloid i,j.
        ! aa          : ai+aj.
        ! Aij         : lubrication correction coefficient for i,j.
        ! Bij         : repulsive force coefficient for i,j.
        !
        ! x_image     : position of closest image of colloid j to i.
        ! v_image     : velocity of closest image of colloid j to i.
        !----------------------------------------------------
        
        REAL(MK)                                :: ai, aj, aa
        REAL(MK)                                :: Aij
        REAL(MK), DIMENSION(3)                  :: vi_old, vj_old
        REAL(MK), DIMENSION(3)                  :: vi_new, vj_new
        
        REAL(MK), DIMENSION(3)                  :: x_image, v_image 
        REAL(MK), DIMENSION(3)                  :: rij, eij, vij
        REAL(MK)                                :: r, h
          
        INTEGER                                 :: i, j
        
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------

        stat_info     = 0
        stat_info_sub = 0
        
        IF ( num_sweep < 1 ) THEN
           PRINT *, __FILE__,__LINE__,&
                "number of sweeps is wrong!"
           stat_info = -1
           GOTO 9999
        END IF

        !----------------------------------------------------
        ! set up physics parameters.
        !----------------------------------------------------
        
        dim   = this%num_dim
        num   = this%num_colloid
        
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
        
        !----------------------------------------------------
        ! set up parameters for implicit sweeps.
        !----------------------------------------------------
        
        dt_sweep  = dt / num_sweep
        
        !----------------------------------------------------
        ! For now, we assume that each radius is the same.
        !----------------------------------------------------
        
        ai = this%radius(1,1)
        aj = ai
        aa = ai+aj
        
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
                    
                 END DO ! j = i + 1, num
                 
              END DO ! i = 1, num - 1
              
           END DO ! sweep = 1, num_sweep
           
        END IF ! cc_lub_type > mcf_cc_lub_type_no
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE colloid_compute_interaction_implicit_velocity_pair_sweep
      
      
      
