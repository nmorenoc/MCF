      SUBROUTINE particles_compute_aeval(this,num,stat_info)
        !----------------------------------------------------
        ! Subroutine :  particles_compute_aeval
        !----------------------------------------------------
        !
        ! Purpose    :  Compute the accleration of egenvalues
        !               of particles, in case of egenvector
        !               dynamics for Non-Newtonian Oldroyd-B
        !               model.
        !
        ! Reference  :  Vazquez-Quesada, Ellero, Espanol
        !               Phyical Review E 79. 056707, 2009.
        !
        ! Remark     :
        !
        ! Revision   :  V0.1  25.08.2010, The number of polymers 
        !               per unit of volume this%n_p is substituted
        !               by the number of polymers per particle,
        !               computed as this%n_p / this%rho(i). The
        !               n_p variable is now real instead of integer.
        !               Deleted an unnecesary line (aeval = 0).
        !               (Adolfo)
        !
        !               V0.1  04.08.2009, original version.
        !
        !----------------------------------------------------
        ! Author    : Xin Bian
        ! Contact   : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !
        ! this           : an object of Particles Class.
        ! num            : number of particles needed to be updated,
        !                   i.e. first num particles in this%x 
        !                   are operated.
        ! lambda         : coefficient required.
        ! stat_info      : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        LOGICAL                                 :: Brownian
        INTEGER                                 :: dim
        REAL(MK)                                :: dt
        REAL(MK)                                :: tau
        REAL(MK)                                :: n_p
        INTEGER                                 :: i,j
        REAL(MK), DIMENSION(:), pointer         :: r_num
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        !----------------------------------------------------
        ! Calculation only for real particles.
        !----------------------------------------------------
        
        IF( num > this%num_part_real) THEN
           PRINT *, "particles_compute_aeval : ", &
                "num > num_part_real, wrong !"
           stat_info = -1
           GOTO 9999
        END IF
        
        Brownian = control_get_Brownian(this%ctrl,stat_info_sub)
        dim = physics_get_num_dim(this%phys,stat_info_sub)
        dt  = physics_get_dt(this%phys,stat_info_sub)
        tau = physics_get_tau(this%phys,stat_info_sub)
        n_p = physics_get_n_p(this%phys,stat_info_sub)
        
        IF ( ABS(tau) < mcf_machine_zero ) THEN
           PRINT *, "particles_compute_aeval : ", &
                "tau shouldn't be zero !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(this%aeval)) THEN
           DEALLOCATE(this%aeval)
        ENd IF
        
        ALLOCATE(this%aeval(dim,num))

!!$        this%aeval(:,:) = 0.0_MK
        
        !------------------------------------------
        ! Compute accleration of egenvalues using 
        ! egenvector dynamics velocity gradient 
        ! tensor and egenvalues themselves.
        !------------------------------------------

        DO j = 1, num ! particle index
           
           DO i =1, dim ! dimension
              
              this%aeval(i,j) = &
                   2.0_MK* this%eval(i,j) * &
                   this%evgt(i+dim*(i-1),j) + &
                   (1.0_MK - this%eval(i,j)) / tau 
                 
           END DO
           
        END DO

        IF ( Brownian ) THEN

           ALLOCATE(r_num(dim))
           
           DO j = 1, num ! particle index

              DO i = 1, dim

                 r_num(i) = random_random(this%random,stat_info_sub)

              ENDDO

              this%aeval(:,j) = this%aeval(:,j) + &
                      2.0_MK/(tau * n_p / this%rho(j)) + &
                      SQRT(4.0_MK*this%eval(:,j) / &
                      (tau * n_p / this%rho(j))) * &
                      r_num(:) / SQRT(dt)
              
           ENDDO                
           
           DEALLOCATE(r_num)

        END IF ! Brownian
        
9999    CONTINUE      
        
        RETURN
        
      END SUBROUTINE particles_compute_aeval
      
      
