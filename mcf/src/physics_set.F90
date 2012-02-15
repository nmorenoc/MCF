!--------------------------------------------------
!  Subroutine   :  physics_set_*
!--------------------------------------------------
!
!  Purpose      : Set routines of Class Physics.
!
!  Reference    :
!
!  Remark       :
!
!  Revisions    : V0.1 01.03.2009, original version.
!
!--------------------------------------------------
! Author       : Xin Bian
! Contact      : xin.bian@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!--------------------------------------------------

      SUBROUTINE physics_set_num_species(this,d_num_species,stat_info)
        
        TYPE(Physics), INTENT(INOUT)            :: this
        INTEGER, INTENT(IN)                     :: d_num_species
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        this%num_species = d_num_species 
        
        RETURN
        
      END SUBROUTINE physics_set_num_species
      
      
      SUBROUTINE physics_set_num_dim(this,d_num_dim,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_num_dim
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: stat_info_sub
        
        stat_info      = 0
        stat_info_sub = 0
        
        IF( d_num_dim < 2 .OR. d_num_dim > 3 ) THEN
           PRINT *, "physics_set_num_dim : ", &
                "Dimension is not supported !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( d_num_dim /= this%num_dim) THEN
           
           this%num_dim = d_num_dim        
           
           IF(ASSOCIATED(this%min_phys)) THEN
              DEALLOCATE(this%min_phys)
           END IF
           ALLOCATE(this%min_phys(d_num_dim))
           
           IF(ASSOCIATED(this%max_phys)) THEN
              DEALLOCATE(this%max_phys)
           END IF
           ALLOCATE(this%max_phys(d_num_dim))

           IF(ASSOCIATED(this%min_phys_t)) THEN
              DEALLOCATE(this%min_phys_t)
           END IF
           ALLOCATE(this%min_phys_t(d_num_dim))
           
           IF(ASSOCIATED(this%max_phys_t)) THEN
              DEALLOCATE(this%max_phys_t)
           END IF
           ALLOCATE(this%max_phys_t(d_num_dim))
           
           IF(ASSOCIATED(this%num_part_dim)) THEN
              DEALLOCATE(this%num_part_dim)
           END IF
           ALLOCATE(this%num_part_dim(d_num_dim))
           
           IF(ASSOCIATED(this%num_part_dim_t)) THEN
              DEALLOCATE(this%num_part_dim_t)
           END IF
           ALLOCATE(this%num_part_dim_t(d_num_dim))
           
           IF(ASSOCIATED(this%dx)) THEN
              DEALLOCATE(this%dx)
           END IF
           ALLOCATE(this%dx(d_num_dim))
           
           IF(ASSOCIATED(this%body_force)) THEN
              DEALLOCATE(this%body_force)
           END IF
           ALLOCATE(this%body_force(d_num_dim))
           
           IF(ASSOCIATED(this%body_force_d)) THEN
              DEALLOCATE(this%body_force_d)
           END IF
           ALLOCATE(this%body_force_d(d_num_dim))
           
           IF(ASSOCIATED(this%bcdef)) THEN
              DEALLOCATE(this%bcdef)
           END IF
           ALLOCATE(this%bcdef(2*d_num_dim))
           
           IF(ASSOCIATED(this%eval)) THEN
              DEALLOCATE(this%eval)
           END IF
           ALLOCATE(this%eval(d_num_dim))
           
           IF(ASSOCIATED(this%evec)) THEN
              DEALLOCATE(this%evec)
           END IF
           ALLOCATE(this%evec(d_num_dim,d_num_dim))
           
           
        END IF
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE physics_set_num_dim
      
      
      SUBROUTINE physics_set_min_phys(this,d_min_phys,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), DIMENSION(:)          :: d_min_phys
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        dim = SIZE(d_min_phys)
        IF (dim /= this%num_dim ) THEN
           PRINT *, "physics_set_min_phys : ", "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%min_phys(1:dim) = d_min_phys(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE physics_set_min_phys
      
      
      SUBROUTINE physics_set_max_phys(this,d_max_phys,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), DIMENSION(:)          :: d_max_phys
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        dim = SIZE(d_max_phys)
        IF (dim /= this%num_dim ) THEN
           PRINT *, "physics_set_max_phys : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999           
        END IF
        
        this%max_phys(1:dim) = d_max_phys(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE physics_set_max_phys
      

      SUBROUTINE physics_set_min_phys_t(this,d_min_phys_t,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), DIMENSION(:)          :: d_min_phys_t
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        dim = SIZE(d_min_phys_t)
        IF (dim /= this%num_dim ) THEN
           PRINT *, "physics_set_min_phys_t : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
9999    CONTINUE
        
        this%min_phys_t(1:dim) = d_min_phys_t(1:dim)
        
        RETURN
        
      END SUBROUTINE physics_set_min_phys_t
      
      
      SUBROUTINE physics_set_max_phys_t(this,d_max_phys_t,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), DIMENSION(:)          :: d_max_phys_t
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        dim = SIZE(d_max_phys_t)
        IF (dim /= this%num_dim ) THEN
           PRINT *, "physics_set_max_phys_t : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%max_phys_t(1:dim) = d_max_phys_t(1:dim)
        
9999    CONTINUE
        RETURN
        
      END SUBROUTINE physics_set_max_phys_t

      
      SUBROUTINE physics_set_lattice_type(this,d_lattice_type,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_lattice_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%lattice_type = d_lattice_type
        
        RETURN
        
      END SUBROUTINE physics_set_lattice_type
      
      
      SUBROUTINE physics_set_num_part_dim(this,d_num_part_dim,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, DIMENSION(:)           :: d_num_part_dim
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_num_part_dim)
        
        IF (dim /= this%num_dim ) THEN
           PRINT *, "physics_set_num_part_dim : ", "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%num_part_dim(1:dim) = d_num_part_dim(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE physics_set_num_part_dim
      
      
      SUBROUTINE physics_set_num_part_dim_t(this,d_num_part_dim_t,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, DIMENSION(:)           :: d_num_part_dim_t
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_num_part_dim_t)
        
        IF (dim /= this%num_dim ) THEN
           PRINT *, "physics_set_num_part_dim_t : ",&
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%num_part_dim_t(1:dim) = d_num_part_dim_t(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE physics_set_num_part_dim_t
      
      
      SUBROUTINE physics_set_num_part_tot(this,d_num_part_tot,stat_info)
       
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_num_part_tot
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        this%num_part_tot = d_num_part_tot
        
        RETURN
        
      END SUBROUTINE physics_set_num_part_tot
      
      
      SUBROUTINE physics_set_dx(this,d_dx,stat_info)
        
        TYPE(Physics), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: d_dx
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        
        stat_info = 0
        
        dim = SIZE(d_dx)
        
        IF (dim /= this%num_dim ) THEN
           PRINT *, "physics_set_dx : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999           
        END IF
        
        this%dx(1:dim) = d_dx(1:dim)      
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE physics_set_dx
      
      
      SUBROUTINE physics_set_cut_off(this,d_cut_off,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_cut_off
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%cut_off = d_cut_off
        
        RETURN
        
      END SUBROUTINE physics_set_cut_off
      
      
      SUBROUTINE physics_set_h(this,d_h,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_h
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%h = d_h
        
        RETURN
        
      END SUBROUTINE physics_set_h
      
      
      SUBROUTINE physics_set_dt_c(this,d_dt,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_dt
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%dt_c = d_dt
        
        RETURN
        
      END SUBROUTINE physics_set_dt_c
      
      
      SUBROUTINE physics_set_dt_nu(this,d_dt,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_dt
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%dt_nu = d_dt
        
        RETURN
        
      END SUBROUTINE physics_set_dt_nu
      
      
      SUBROUTINE physics_set_fa_max(this,d_fa,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_fa
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%fa_max = d_fa
        
        RETURN
        
      END SUBROUTINE physics_set_fa_max
      
      
      SUBROUTINE physics_set_dt_f(this,d_dt,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_dt
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%dt_f = d_dt
        
        RETURN
        
      END SUBROUTINE physics_set_dt_f


      SUBROUTINE physics_set_dt(this,d_dt,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_dt
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%dt = d_dt
        
        RETURN
        
      END SUBROUTINE physics_set_dt
      
      
      SUBROUTINE physics_set_step_start(this,d_step,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_step
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%step_start = d_step
        
        RETURN
        
      END SUBROUTINE physics_set_step_start
      
      
      SUBROUTINE physics_set_step_current(this,d_step,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_step
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%step_current = d_step
        
        RETURN
        
      END SUBROUTINE physics_set_step_current
      

      SUBROUTINE physics_set_step_end(this,d_step_end,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_step_end
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%step_end = d_step_end 
        
        RETURN
        
      END SUBROUTINE physics_set_step_end
      
      
      SUBROUTINE physics_set_time_start(this,d_time,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_time
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%time_start = d_time
        
        RETURN
        
      END SUBROUTINE physics_set_time_start
      
      
      SUBROUTINE physics_set_time_current(this,d_time,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_time
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%time_current = d_time
        
        RETURN
        
      END SUBROUTINE physics_set_time_current

      
      SUBROUTINE physics_set_time_end(this,d_time,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_time
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%time_end = d_time
        
        RETURN
        
      END SUBROUTINE physics_set_time_end
      
      
      SUBROUTINE physics_set_rho(this,d_rho,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_rho                    
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%rho = d_rho
        
        RETURN
        
      END SUBROUTINE physics_set_rho


      SUBROUTINE physics_set_eta(this,d_eta,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_eta
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%eta =  d_eta 
        
        RETURN
        
      END SUBROUTINE physics_set_eta
      

      SUBROUTINE physics_set_eta_coef(this,d_eta_coef,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_eta_coef
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%eta_coef =  d_eta_coef
        
        RETURN
        
      END SUBROUTINE physics_set_eta_coef

      
      SUBROUTINE physics_set_ksai(this,d_ksai,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_ksai
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%ksai = d_ksai
        
        RETURN
        
      END SUBROUTINE physics_set_ksai
      
      
      SUBROUTINE physics_set_kt(this,d_kt,stat_info)
        
        TYPE(Physics), INTENT(INOUT)     :: this
        REAL(MK), INTENT(IN)             :: d_kt
        INTEGER, INTENT(OUT)             :: stat_info
        
        stat_info = 0
        
        this%kt = d_kt
        
        RETURN
        
      END SUBROUTINE physics_set_kt
      
      
      SUBROUTINE physics_set_c(this,d_c,stat_info)
        
        TYPE(Physics), INTENT(INOUT)     :: this
        REAL(MK), INTENT(IN)             :: d_c
        INTEGER, INTENT(OUT)             :: stat_info
        
        stat_info = 0
        
        this%c = d_c
        
        RETURN
        
      END SUBROUTINE physics_set_c
      
      
      SUBROUTINE physics_set_rho_ref(this,d_rho_ref,stat_info)
        
        TYPE(Physics), INTENT(INOUT)     :: this
        REAL(MK), INTENT(IN)             :: d_rho_ref
        INTEGER, INTENT(OUT)             :: stat_info
        
        stat_info = 0
        
        this%rho_ref = d_rho_ref
        
        RETURN
        
      END SUBROUTINE physics_set_rho_ref
      
      
      SUBROUTINE physics_set_gamma(this,d_gamma,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_gamma
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%gamma = d_gamma
        
        RETURN
        
      END SUBROUTINE physics_set_gamma
      
      
      SUBROUTINE physics_set_relax_type(this,d_relax_type,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_relax_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%relax_type = d_relax_type
        
        RETURN
        
      END SUBROUTINE physics_set_relax_type
      

      SUBROUTINE physics_set_dt_relax(this,d_dt,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_dt
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%dt_relax = d_dt
        
        RETURN
        
      END SUBROUTINE physics_set_dt_relax


      SUBROUTINE physics_set_step_relax(this,d_step_relax,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_step_relax
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%step_relax = d_step_relax
        
        RETURN
        
      END SUBROUTINE physics_set_step_relax
      

      SUBROUTINE physics_set_time_relax(this,d_time_relax,stat_info)
        
        TYPE(Physics), INTENT(INOUT)     :: this
        REAL(MK), INTENT(IN)             :: d_time_relax
        INTEGER, INTENT(OUT)             :: stat_info
        
        stat_info = 0
        
        this%time_relax = d_time_relax
        
        RETURN
        
      END SUBROUTINE physics_set_time_relax
      

      SUBROUTINE physics_set_disorder_level(this,d_disorder_level,stat_info)
        
        TYPE(Physics), INTENT(INOUT)     :: this
        REAL(MK), INTENT(IN)             :: d_disorder_level
        INTEGER, INTENT(OUT)             :: stat_info
        
        stat_info = 0
        
        this%disorder_level = d_disorder_level
        
        RETURN
        
      END SUBROUTINE physics_set_disorder_level

      
      SUBROUTINE physics_set_kt_relax(this,d_kt_relax,stat_info)
        
        TYPE(Physics), INTENT(INOUT)     :: this
        REAL(MK), INTENT(IN)             :: d_kt_relax
        INTEGER, INTENT(OUT)             :: stat_info
        
        stat_info = 0
        
        this%kt_relax = d_kt_relax
        
        RETURN
        
      END SUBROUTINE physics_set_kt_relax
      
      
      SUBROUTINE physics_set_c_relax(this,d_c_relax,stat_info)
         
        TYPE(Physics), INTENT(INOUT)     :: this
        REAL(MK), INTENT(IN)             :: d_c_relax
        INTEGER, INTENT(OUT)             :: stat_info
        
        stat_info = 0
        
        this%c_relax = d_c_relax
        
        RETURN
        
      END SUBROUTINE physics_set_c_relax


      SUBROUTINE physics_set_tau(this,d_tau,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_tau
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%tau =  d_tau
        
        RETURN
        
      END SUBROUTINE physics_set_tau
      
      SUBROUTINE physics_set_tau_sm(this,d_tau_sm,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_tau_sm
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%tau_sm =  d_tau_sm
        
        RETURN
        
      END SUBROUTINE physics_set_tau_sm

      SUBROUTINE physics_set_k_sm(this,d_k_sm,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_k_sm
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%k_sm =  d_k_sm
        
        RETURN
        
      END SUBROUTINE physics_set_k_sm

      SUBROUTINE physics_set_n_p(this,d_n_p,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_n_p
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%n_p =  d_n_p
        
        RETURN
        
      END SUBROUTINE physics_set_n_p
      
        
      SUBROUTINE physics_set_kt_p(this,d_kt_p,stat_info)
        
        TYPE(Physics), INTENT(INOUT)     :: this
        REAL(MK), INTENT(IN)             :: d_kt_p
        INTEGER, INTENT(OUT)             :: stat_info
        
        stat_info = 0
        
        this%kt_p = d_kt_p
        
        RETURN
        
      END SUBROUTINE physics_set_kt_p
      
      
      SUBROUTINE physics_set_eigen_dynamics(this,d_eigen_dynamics,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        LOGICAL, INTENT(IN)             :: d_eigen_dynamics
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%eigen_dynamics =  d_eigen_dynamics
        
        RETURN
        
      END SUBROUTINE physics_set_eigen_dynamics
   

      SUBROUTINE physics_set_eval(this,d_eval,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), DIMENSION(:)          :: d_eval
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        
        stat_info = 0
        
        dim = SIZE(d_eval)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "physics_set_eval : ", "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%eval(1:dim) = d_eval(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  physics_set_eval
      
      
      SUBROUTINE physics_set_eval_tolerance(this,d_tol,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK)                        :: d_tol
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        this%eval_tolerance = d_tol
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  physics_set_eval_tolerance
      

      SUBROUTINE physics_set_evec(this,d_evec,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), DIMENSION(:,:)        :: d_evec
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER, DIMENSION(2)           :: dim
        
        
        stat_info = 0
        
        dim(1) = SIZE(d_evec,1)
        dim(2) = SIZE(d_evec,2)
        
        IF( dim(1) /= this%num_dim .OR. &
             dim(2) /= this%num_dim) THEN
           PRINT *, "physics_set_evec : ", "Wrong dimension !"
           stat_info = -1
           GOTO 9999           
        END IF
        
        this%evec(1:dim(1),1:dim(2)) = &
             d_evec(1:dim(1),1:dim(2))
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  physics_set_evec
      

      SUBROUTINE physics_set_evec_normalize(this,d_norm,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        LOGICAL                         :: d_norm
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        this%evec_normalize = d_norm
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  physics_set_evec_normalize

      
      SUBROUTINE physics_set_evec_tolerance(this,d_tol,stat_info)
         
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK)                        :: d_tol
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        this%evec_tolerance = d_tol
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  physics_set_evec_tolerance

      
      SUBROUTINE physics_set_body_force_type(this,d_body_force_type,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_body_force_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%body_force_type =  d_body_force_type
        
        RETURN
        
      END SUBROUTINE physics_set_body_force_type
      
      
      SUBROUTINE physics_set_body_force(this,d_body_force,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), DIMENSION(:)          :: d_body_force
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_body_force)
        
        IF ( dim /= this%num_dim) THEN
           PRINT *, "physics_set_body_foce : ", "Wrong Dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%body_force(1:dim) = d_body_force(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  physics_set_body_force
      
      
      SUBROUTINE physics_set_body_force_d(this,d_body_force_d,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), DIMENSION(:)          :: d_body_force_d
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_body_force_d)
        
        IF ( dim /= this%num_dim) THEN
           PRINT *, "physics_set_body_foce_d : ", "Wrong Dimension !"
           stat_info = -1
           GOTO 9999           
        END IF

        this%body_force_d(1:dim) = d_body_force_d(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  physics_set_body_force_d
      
      
      SUBROUTINE physics_set_flow_direction(this,d_flow_direction,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_flow_direction
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%flow_direction = d_flow_direction
        
        RETURN
        
      END SUBROUTINE physics_set_flow_direction


      SUBROUTINE physics_set_flow_width(this,d_flow_width,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_flow_width
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%flow_width = d_flow_width
        
        RETURN
        
      END SUBROUTINE physics_set_flow_width
      
      
      SUBROUTINE physics_set_flow_v(this,d_flow_v,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_flow_v
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%flow_v = d_flow_v
        
        RETURN
        
      END SUBROUTINE physics_set_flow_v
      
  
      SUBROUTINE physics_set_flow_adjust_freq(this,d_flow_adjust_freq,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_flow_adjust_freq
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%flow_adjust_freq = d_flow_adjust_freq
        
        RETURN
        
      END SUBROUTINE physics_set_flow_adjust_freq


      SUBROUTINE physics_set_num_colloid(this,d_num_colloid,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_num_colloid
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%num_colloid =  d_num_colloid
        
        RETURN
        
      END SUBROUTINE physics_set_num_colloid
      
      
      SUBROUTINE physics_set_colloid(this,d_colloid,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        TYPE(Colloid), TARGET           :: d_colloid
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_dim
        
        stat_info     = 0
        stat_info_sub = 0
        
        num_dim = colloid_get_num_dim(d_colloid,stat_info_sub)
        
        IF(num_dim /= this%num_dim) THEN
           PRINT *, "physics_set_colloid : ",&
                "Colloids dimension doesn't match with physics !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%colloids => d_colloid
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE physics_set_colloid
      
      
      SUBROUTINE physics_set_bcdef(this,d_bcdef,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        INTEGER, DIMENSION(:)           :: d_bcdef
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_bcdef)
        
        IF( dim /= 2*this%num_dim) THEN
           PRINT *, "physics_set_bcdef : ", "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%bcdef(1:dim) = d_bcdef(1:dim)

9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  physics_set_bcdef
      

      SUBROUTINE physics_set_boundary(this,d_boundary,stat_info)
        
        TYPE(Physics), INTENT(INOUT)    :: this
        TYPE(Boundary), TARGET          :: d_boundary
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_dim
        
        stat_info     = 0
        stat_info_sub = 0
        
        num_dim = boundary_get_num_dim(d_boundary,stat_info_sub)
        
        IF(num_dim /= this%num_dim) THEN
           PRINT *, "physics_set_boundary : ",&
                "Boundarys dimension doesn't match with physics !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%boundary => d_boundary
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE physics_set_boundary
      
