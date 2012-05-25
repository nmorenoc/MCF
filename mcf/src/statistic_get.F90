!------------------------------------------------------------
! All the public "get" subroutines of Class statistic,
! which return member variables of statistic.
!
! Reference   :
!
! Remark      :
!
! Revisions   :  V0.1 03.03.2009, original version.
!
!------------------------------------------------------------
! Author       : Xin Bian
! Contact      : xin.bian@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!------------------------------------------------------------      

      LOGICAL FUNCTION statistic_get_flow_v_fixed(this,stat_info)

        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        
        statistic_get_flow_v_fixed = this%flow_v_fixed
        
        RETURN
        
      END FUNCTION statistic_get_flow_v_fixed
      
      
      LOGICAL FUNCTION statistic_get_l_p_energy(this,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        
        statistic_get_l_p_energy = this%l_p_energy
        
        RETURN
        
      END FUNCTION statistic_get_l_p_energy
      

      INTEGER FUNCTION statistic_get_num_dim(this,stat_info)

        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        
        statistic_get_num_dim = this%num_dim
        
        RETURN
        
      END FUNCTION  statistic_get_num_dim
      
      
      SUBROUTINE statistic_get_statistic(this, &
           d_k_energy,d_momentum,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        REAL(MK), INTENT(OUT)           :: d_k_energy
        REAL(MK), DIMENSION(:), POINTER :: d_momentum
        INTEGER, INTENT(OUT)            :: stat_info 
        
        
        stat_info = 0
        
        d_k_energy = this%k_energy
        
        IF(ASSOCIATED(d_momentum)) THEN
           DEALLOCATE(d_momentum)
        END IF
        
        
        ALLOCATE(d_momentum(this%num_dim))
        
        d_momentum = this%momentum(1:this%num_dim)
        
        RETURN        
        
      END SUBROUTINE statistic_get_statistic

      
      REAL(MK) FUNCTION statistic_get_kinetic_energy(this,stat_info)

        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        statistic_get_kinetic_energy = this%k_energy
        
        RETURN
        
      END FUNCTION statistic_get_kinetic_energy


      SUBROUTINE statistic_get_momentum(this, &
           d_momentum,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        REAL(MK), DIMENSION(:), POINTER :: d_momentum
        INTEGER, INTENT(OUT)            :: stat_info 
        
        
        stat_info = 0
        
        IF(ASSOCIATED(d_momentum)) THEN
           DEALLOCATE(d_momentum)
        END IF
        
        ALLOCATE(d_momentum(this%num_dim))
        
        d_momentum(:) = this%momentum(1:this%num_dim)
        
        RETURN        
        
      END SUBROUTINE statistic_get_momentum


      SUBROUTINE statistic_get_stress(this, &
           d_stress,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        REAL(MK), DIMENSION(:), POINTER :: d_stress
        INTEGER, INTENT(OUT)            :: stat_info 
        
        INTEGER                         :: dim2

        
        stat_info = 0
        dim2      = this%num_dim**2
        
        IF(ASSOCIATED(d_stress)) THEN
           DEALLOCATE(d_stress)
        END IF
        
        ALLOCATE(d_stress(dim2))
        
        d_stress(1:dim2) = this%stress(1:dim2)
        
        RETURN        
        
      END SUBROUTINE statistic_get_stress
      

      SUBROUTINE statistic_get_stress_p(this, &
           d_stress,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        REAL(MK), DIMENSION(:), POINTER :: d_stress
        INTEGER, INTENT(OUT)            :: stat_info 
        
        INTEGER                         :: dim2

        
        stat_info = 0
        dim2      = this%num_dim**2
        
        IF(ASSOCIATED(d_stress)) THEN
           DEALLOCATE(d_stress)
        END IF
        
        ALLOCATE(d_stress(dim2))
        
        d_stress(1:dim2) = this%stress_p(1:dim2)
        
        RETURN        
        
      END SUBROUTINE statistic_get_stress_p

      
      SUBROUTINE statistic_get_stress_v(this, &
           d_stress,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        REAL(MK), DIMENSION(:), POINTER :: d_stress
        INTEGER, INTENT(OUT)            :: stat_info 
        
        INTEGER                         :: dim2

        
        stat_info = 0
        dim2      = this%num_dim**2
        
        IF(ASSOCIATED(d_stress)) THEN
           DEALLOCATE(d_stress)
        END IF
        
        ALLOCATE(d_stress(dim2))
        
        d_stress(1:dim2) = this%stress_v(1:dim2)
        
        RETURN        
        
      END SUBROUTINE statistic_get_stress_v
      

      SUBROUTINE statistic_get_stress_r(this, &
           d_stress,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        REAL(MK), DIMENSION(:), POINTER :: d_stress
        INTEGER, INTENT(OUT)            :: stat_info 
        
        INTEGER                         :: dim2

        
        stat_info = 0
        dim2      = this%num_dim**2
        
        IF(ASSOCIATED(d_stress)) THEN
           DEALLOCATE(d_stress)
        END IF
        
        ALLOCATE(d_stress(dim2))
        
        d_stress(1:dim2) = this%stress_r(1:dim2)
        
        RETURN        
        
      END SUBROUTINE statistic_get_stress_r

      
      SUBROUTINE statistic_get_disorder(this,d_disorder,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        REAL(MK), DIMENSION(:), POINTER :: d_disorder
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0        
        
        IF(ASSOCIATED(d_disorder)) THEN
           DEALLOCATE(d_disorder)
        END IF        
        
        ALLOCATE(d_disorder(this%num_dim))
        d_disorder  = this%disorder(1:this%num_dim)
        
        RETURN        
        
      END SUBROUTINE statistic_get_disorder
      
      
      REAL(MK) FUNCTION statistic_get_disorder_max(this,stat_info)

        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        REAL(MK)                        :: disorder_max
        INTEGER                         :: i
        
        stat_info = 0
        
        disorder_max = this%disorder(1)

        DO i = 2, this%num_dim
           
           IF ( this%disorder(i) > disorder_max ) THEN
              
              disorder_max = this%disorder(i)
              
           END IF
           
        END DO
        
        statistic_get_disorder_max = disorder_max
        
        RETURN
        
      END FUNCTION statistic_get_disorder_max

 
      REAL(MK) FUNCTION statistic_get_disorder_square_root(this,stat_info)

        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        REAL(MK)                        :: disorder
        INTEGER                         :: i
        
        stat_info = 0
        
        disorder = 0.0_MK

        DO i = 1, this%num_dim
           
           disorder = disorder + this%disorder(i)**2
           
        END DO
        
        statistic_get_disorder_square_root = SQRT(disorder)
        
        RETURN
        
      END FUNCTION statistic_get_disorder_square_root


      REAL(MK) FUNCTION statistic_get_disorder_average(this,stat_info)

        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        REAL(MK)                        :: disorder
        INTEGER                         :: i
        
        stat_info = 0
        
        disorder = 0.0_MK

        DO i = 1, this%num_dim
           
           disorder = disorder + this%disorder(i)
           
        END DO
        
        statistic_get_disorder_average = disorder / this%num_dim
        
        RETURN
        
      END FUNCTION statistic_get_disorder_average
      
      
      SUBROUTINE statistic_get_v_average(this,d_v_aver,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        REAL(MK), DIMENSION(:), POINTER :: d_v_aver
        INTEGER, INTENT(OUT)            :: stat_info 
        
        
        stat_info = 0        
        
        IF(ASSOCIATED(d_v_aver)) THEN
           DEALLOCATE(d_v_aver)
        END IF
        
        ALLOCATE(d_v_aver(this%num_dim))
        
        d_v_aver(1:this%num_dim) = &
             this%v_aver(1:this%num_dim)
        
        RETURN        
        
      END SUBROUTINE statistic_get_v_average
      
      
      REAL(MK) FUNCTION statistic_get_rho_min(this,stat_info)

        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        
        statistic_get_rho_min = this%rho_min
        
        RETURN
        
      END FUNCTION  statistic_get_rho_min

      
      REAL(MK) FUNCTION statistic_get_rho_max(this,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        
        statistic_get_rho_max = this%rho_max
        
        RETURN
        
      END FUNCTION  statistic_get_rho_max


      REAL(MK) FUNCTION statistic_get_p_energy(this,stat_info)

        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        
        statistic_get_p_energy = this%p_energy
        
        RETURN
        
      END FUNCTION statistic_get_p_energy


      LOGICAL FUNCTION statistic_get_colloid_implicit_pair_sweep_adaptive(this,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        
        statistic_get_colloid_implicit_pair_sweep_adaptive = &
             this%colloid_implicit_pair_sweep_adaptive
        
        RETURN        
        
      END FUNCTION statistic_get_colloid_implicit_pair_sweep_adaptive

      
      INTEGER FUNCTION statistic_get_colloid_implicit_pair_num_sweep(this,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        
        statistic_get_colloid_implicit_pair_num_sweep = &
             this%colloid_implicit_pair_num_sweep 
        
        RETURN        
        
      END FUNCTION  statistic_get_colloid_implicit_pair_num_sweep
      
      
      REAL(MK) FUNCTION statistic_get_colloid_implicit_pair_sweep_error(this,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        
        statistic_get_colloid_implicit_pair_sweep_error = &
             this%colloid_implicit_pair_sweep_error 
        
        RETURN        
        
      END FUNCTION statistic_get_colloid_implicit_pair_sweep_error
      
