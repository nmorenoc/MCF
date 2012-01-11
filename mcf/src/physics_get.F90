!--------------------------------------------------
! Subroutine  : physics_get_*
!--------------------------------------------------
!
! Purpose     : Get routines of Class Physics.
!
! Reference   :
!
! Remark      :
!
! Revisions   : V0.1 01.03.2009, original version.
!
!--------------------------------------------------
! Author      : Xin Bian
! Contact     : xin.bian@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!--------------------------------------------------

      INTEGER FUNCTION physics_get_num_species(this,stat_info)
       
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_num_species = this%num_species
        
        RETURN
        
      END FUNCTION physics_get_num_species
      
      
      INTEGER FUNCTION physics_get_num_dim(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_num_dim = this%num_dim
        
        RETURN
        
      END FUNCTION physics_get_num_dim
      
      
      SUBROUTINE physics_get_min_phys(this,d_min_phys,stat_info)
       
        TYPE(Physics), INTENT(IN)       :: this
        REAL(MK), DIMENSION(:), POINTER :: d_min_phys
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_min_phys))THEN
           DEALLOCATE(d_min_phys)
        END IF
        
        ALLOCATE(d_min_phys(this%num_dim))
        d_min_phys(1:this%num_dim) = this%min_phys(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE physics_get_min_phys
      
      
      SUBROUTINE physics_get_max_phys(this,d_max_phys,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        REAL(MK), DIMENSION(:), POINTER :: d_max_phys
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_max_phys))THEN
           DEALLOCATE(d_max_phys)
        END IF
        
        ALLOCATE(d_max_phys(this%num_dim))
        d_max_phys(1:this%num_dim) = this%max_phys(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE physics_get_max_phys
      

      SUBROUTINE physics_get_length(this,d_length,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        REAL(MK), DIMENSION(:), POINTER :: d_length
        INTEGER, INTENT(OUT)            :: stat_info

        
        stat_info = 0
        
        IF(ASSOCIATED(d_length))THEN
           DEALLOCATE(d_length)
        END IF
        
        ALLOCATE(d_length(1:this%num_dim))
        
        d_length(:) = &
             this%max_phys(:) - this%min_phys(:)
        
        RETURN
        
      END SUBROUTINE physics_get_length

      
      SUBROUTINE physics_get_min_phys_t(this,d_min_phys_t,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        REAL(MK), DIMENSION(:), POINTER :: d_min_phys_t
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_min_phys_t))THEN
           DEALLOCATE(d_min_phys_t)
        END IF
        
        ALLOCATE(d_min_phys_t(this%num_dim))
        d_min_phys_t(1:this%num_dim) = this%min_phys_t(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE physics_get_min_phys_t
      
      
      SUBROUTINE physics_get_max_phys_t(this,d_max_phys_t,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        REAL(MK), DIMENSION(:), POINTER :: d_max_phys_t
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_max_phys_t))THEN
           DEALLOCATE(d_max_phys_t)
        END IF
        
        ALLOCATE(d_max_phys_t(this%num_dim))
        d_max_phys_t(1:this%num_dim) = this%max_phys_t(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE physics_get_max_phys_t
      
      
      SUBROUTINE physics_get_length_t(this,d_length_t,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        REAL(MK), DIMENSION(:), POINTER :: d_length_t
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_length_t))THEN
           DEALLOCATE(d_length_t)
        END IF
        
        ALLOCATE(d_length_t(this%num_dim))
        
        d_length_t(1:this%num_dim) = &
             this%max_phys_t(1:this%num_dim) - &
             this%min_phys_t(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE physics_get_length_t
      
      
      INTEGER FUNCTION physics_get_lattice_type(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_lattice_type = this%lattice_type
        
        RETURN
        
      END FUNCTION physics_get_lattice_type
      
      
      SUBROUTINE physics_get_num_part_dim(this,d_num_part_dim,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, DIMENSION(:), POINTER  :: d_num_part_dim
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_num_part_dim))THEN
           DEALLOCATE(d_num_part_dim)
        END IF
        
        ALLOCATE(d_num_part_dim(this%num_dim))
        d_num_part_dim(1:this%num_dim) = this%num_part_dim(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE physics_get_num_part_dim

      
      SUBROUTINE physics_get_num_part_dim_t(this,d_num_part_dim_t,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, DIMENSION(:), POINTER  :: d_num_part_dim_t
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_num_part_dim_t))THEN
           DEALLOCATE(d_num_part_dim_t)
        END IF
        
        ALLOCATE(d_num_part_dim_t(this%num_dim))
        d_num_part_dim_t(1:this%num_dim) = &
             this%num_part_dim_t(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE physics_get_num_part_dim_t
      
      
      INTEGER FUNCTION physics_get_num_part_tot(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_num_part_tot = this%num_part_tot
        
        RETURN
        
      END FUNCTION physics_get_num_part_tot
      
      
      SUBROUTINE physics_get_dx(this,d_dx,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        REAL(MK), DIMENSION(:), POINTER :: d_dx
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        IF(ASSOCIATED(d_dx))THEN
           DEALLOCATE(d_dx)
        END IF
        
        ALLOCATE(d_dx(this%num_dim))
        d_dx(1:this%num_dim) = this%dx(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE physics_get_dx
      
      
      REAL(MK) FUNCTION physics_get_cut_off(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_cut_off = this%cut_off
        
        RETURN
        
      END FUNCTION physics_get_cut_off
      
      
      REAL(MK) FUNCTION physics_get_h(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_h = this%h
        
        RETURN
        
      END FUNCTION physics_get_h
      
      
       REAL(MK) FUNCTION physics_get_dt_c(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_dt_c = this%dt_c
        
        RETURN
        
      END FUNCTION physics_get_dt_c

      
      REAL(MK) FUNCTION physics_get_dt_nu(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_dt_nu = this%dt_nu
        
        RETURN
        
      END FUNCTION physics_get_dt_nu
      
      
      REAL(MK) FUNCTION physics_get_fa_max(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_fa_max = this%fa_max
        
        RETURN
        
      END FUNCTION physics_get_fa_max

      
      REAL(MK) FUNCTION physics_get_dt_f(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_dt_f = this%dt_f
        
        RETURN
        
      END FUNCTION physics_get_dt_f


      REAL(MK) FUNCTION physics_get_dt(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_dt = this%dt
        
        RETURN
        
      END FUNCTION physics_get_dt
      
      
      INTEGER FUNCTION physics_get_step_start(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_step_start = this%step_start
        
        RETURN
        
      END FUNCTION physics_get_step_start
      

      INTEGER FUNCTION physics_get_step_current(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_step_current = this%step_current
        
        RETURN
        
      END FUNCTION physics_get_step_current


      INTEGER FUNCTION physics_get_step_end(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_step_end = this%step_end
        
        RETURN
        
      END FUNCTION physics_get_step_end
      
      
      REAL(MK) FUNCTION physics_get_time_start(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_time_start = this%time_start
        
        RETURN
        
      END FUNCTION physics_get_time_start
      

      REAL(MK) FUNCTION physics_get_time_current(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_time_current = this%time_current
        
        RETURN
        
      END FUNCTION physics_get_time_current
      

      REAL(MK) FUNCTION physics_get_time_end(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_time_end = this%time_end
        
        RETURN
        
      END FUNCTION physics_get_time_end
      
      
      REAL(MK) FUNCTION physics_get_rho(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_rho = this%rho
        
        RETURN
        
      END FUNCTION physics_get_rho
      

      REAL(MK) FUNCTION physics_get_eta(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_eta = this%eta
        
        RETURN
        
      END FUNCTION physics_get_eta

      
      REAL(MK) FUNCTION physics_get_eta_coef(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_eta_coef = this%eta_coef
        
        RETURN
        
      END FUNCTION physics_get_eta_coef

       
      REAL(MK) FUNCTION physics_get_ksai(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_ksai = this%ksai
        
        RETURN
        
      END FUNCTION physics_get_ksai
      

      REAL(MK) FUNCTION physics_get_kt(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_kt = this%kt
        
        RETURN
        
      END FUNCTION physics_get_kt

        
      REAL(MK) FUNCTION physics_get_c(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_c = this%c
        
        RETURN
        
      END FUNCTION physics_get_c

                     
      REAL(MK) FUNCTION physics_get_rho_ref(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_rho_ref = this%rho_ref
        
        RETURN
        
      END FUNCTION physics_get_rho_ref
      
      
      REAL(MK) FUNCTION physics_get_gamma(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_gamma = this%gamma
        
        RETURN
        
      END FUNCTION physics_get_gamma
      
      
      INTEGER FUNCTION physics_get_relax_type(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_relax_type = this%relax_type
        
        RETURN
        
      END FUNCTION physics_get_relax_type
      
      
      REAL(MK) FUNCTION physics_get_dt_relax(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_dt_relax = this%dt_relax
        
        RETURN
        
      END FUNCTION physics_get_dt_relax
      
      
      INTEGER FUNCTION physics_get_step_relax(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_step_relax = this%step_relax
        
        RETURN
        
      END FUNCTION physics_get_step_relax
      
      
      REAL(MK) FUNCTION physics_get_time_relax(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_time_relax = this%time_relax
        
        RETURN
        
      END FUNCTION physics_get_time_relax


      REAL(MK) FUNCTION physics_get_disorder_level(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_disorder_level = this%disorder_level
        
        RETURN
        
      END FUNCTION physics_get_disorder_level
      
      
      REAL(MK) FUNCTION physics_get_kt_relax(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_kt_relax = this%kt_relax
        
        RETURN
        
      END FUNCTION physics_get_kt_relax


       REAL(MK) FUNCTION physics_get_c_relax(this,stat_info)
         
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_c_relax = this%c_relax
        
        RETURN
        
      END FUNCTION physics_get_c_relax

      
      REAL(MK) FUNCTION physics_get_tau(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_tau = this%tau
        
        RETURN
        
      END FUNCTION physics_get_tau
      
      REAL(MK) FUNCTION physics_get_tau_sm(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_tau_sm = this%tau_sm
        
        RETURN
        
      END FUNCTION physics_get_tau_sm

      REAL(MK) FUNCTION physics_get_k_sm(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_k_sm = this%k_sm
        
        RETURN
        
      END FUNCTION physics_get_k_sm
   
      REAL(MK) FUNCTION physics_get_n_p(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_n_p = this%n_p
        
        RETURN
        
      END FUNCTION physics_get_n_p
      
      
      REAL(MK) FUNCTION physics_get_kt_p(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_kt_p = this%kt_p
        
        RETURN
        
      END FUNCTION physics_get_kt_p

      
      LOGICAL FUNCTION physics_get_eigen_dynamics(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_eigen_dynamics = this%eigen_dynamics
        
        RETURN
        
      END FUNCTION physics_get_eigen_dynamics
      
   
      SUBROUTINE physics_get_eval(this,d_eval,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        REAL(MK), DIMENSION(:), POINTER :: d_eval
        INTEGER, INTENT(OUT)            :: stat_info
          
        stat_info = 0
        
        IF(ASSOCIATED(d_eval)) THEN
           DEALLOCATE(d_eval)
        END IF
        
        ALLOCATE(d_eval(this%num_dim))
        
        d_eval(1:this%num_dim) = this%eval(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE  physics_get_eval
      

      REAL(MK) FUNCTION physics_get_eval_tolerance(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_eval_tolerance = this%eval_tolerance
        
        RETURN
        
      END FUNCTION physics_get_eval_tolerance

      
      SUBROUTINE physics_get_evec(this,d_evec,stat_info)
        
        TYPE(Physics), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_evec
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_evec)) THEN
           DEALLOCATE(d_evec)
        END IF
        
        ALLOCATE(d_evec(this%num_dim,this%num_dim))
        
        d_evec(1:this%num_dim,1:this%num_dim) = &
             this%evec(1:this%num_dim,1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE  physics_get_evec
      

      LOGICAL FUNCTION physics_get_evec_normalize(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_evec_normalize = this%evec_normalize
        
        RETURN
        
      END FUNCTION physics_get_evec_normalize
      

      REAL(MK) FUNCTION physics_get_evec_tolerance(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_evec_tolerance = this%evec_tolerance
        
        RETURN
        
      END FUNCTION physics_get_evec_tolerance
      
      
      INTEGER FUNCTION physics_get_body_force_type(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_body_force_type = this%body_force_type
        
        RETURN
        
      END FUNCTION physics_get_body_force_type


      SUBROUTINE physics_get_body_force(this,d_body_force,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        REAL(MK), DIMENSION(:), POINTER :: d_body_force
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_body_force)) THEN
           DEALLOCATE(d_body_force)
        END IF
        
        ALLOCATE(d_body_force(this%num_dim))
        
        d_body_force(1:this%num_dim) = this%body_force(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE  physics_get_body_force


      SUBROUTINE physics_get_body_force_d(this,d_body_force_d,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        REAL(MK), DIMENSION(:), POINTER :: d_body_force_d
        INTEGER, INTENT(OUT)            :: stat_info
          
        stat_info = 0
        
        IF(ASSOCIATED(d_body_force_d)) THEN
           DEALLOCATE(d_body_force_d)
        END IF
        
        ALLOCATE(d_body_force_d(this%num_dim))
        
        d_body_force_d(1:this%num_dim) = this%body_force_d(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE physics_get_body_force_d
      
      
      INTEGER FUNCTION physics_get_flow_direction(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_flow_direction = this%flow_direction
        
        RETURN
        
      END FUNCTION physics_get_flow_direction

      
      REAL(MK) FUNCTION physics_get_flow_width(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_flow_width = this%flow_width
        
        RETURN
        
      END FUNCTION physics_get_flow_width
      
  
      REAL(MK) FUNCTION physics_get_flow_v(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_flow_v = this%flow_v
        
        RETURN
        
      END FUNCTION physics_get_flow_v
      
        
      INTEGER FUNCTION physics_get_flow_adjust_freq(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_flow_adjust_freq = this%flow_adjust_freq
        
        RETURN
        
      END FUNCTION physics_get_flow_adjust_freq
      

      INTEGER FUNCTION physics_get_num_colloid(this,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        physics_get_num_colloid = this%num_colloid
        
        RETURN
        
      END FUNCTION physics_get_num_colloid
      
      
      SUBROUTINE physics_get_colloid(this,d_colloid,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        TYPE(Colloid), POINTER          :: d_colloid
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        d_colloid => this%colloids
        
        RETURN
        
      END SUBROUTINE physics_get_colloid


      SUBROUTINE physics_get_bcdef(this,d_bcdef,stat_info)
        
        TYPE(Physics), INTENT(IN)       :: this
        INTEGER, DIMENSION(:), POINTER  :: d_bcdef
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_bcdef)) THEN
           DEALLOCATE(d_bcdef)
        END IF
        
        ALLOCATE(d_bcdef(2*this%num_dim))
        
        d_bcdef(1:(2*this%num_dim)) = this%bcdef(1:(2*this%num_dim))
        
        RETURN
        
      END SUBROUTINE  physics_get_bcdef

        
      SUBROUTINE physics_get_boundary(this,d_boundary,stat_info)
    
        TYPE(Physics), INTENT(IN)       :: this
        TYPE(Boundary), POINTER         :: d_boundary
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        d_boundary => this%boundary
        
        RETURN
        
      END SUBROUTINE physics_get_boundary
