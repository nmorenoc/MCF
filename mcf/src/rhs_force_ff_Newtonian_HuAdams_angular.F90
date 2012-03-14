!----------------------------------------------------------------
! This contains the routines which calculate the force and 
! return two accelerations between two Newtonian particles, 
! accelerations from i.e. conservative(pressure) force 
! dissipative(viscous) force and random(thermal noise) force.
! Using Hu & Adams (2006 Physics of fluids) formulation.
!----------------------------------------------------------------
#ifdef __PARTICLES_FORCE_SEPARATE
      SUBROUTINE rhs_force_ff_Newtonian_HuAdams_angular(this,&
           xi,xj,dij,vi,vj,numi,numj,pi,pj,&
           mi,mj,w,gradw,fi,fj,fpi,fpj,fvi,fvj,fri,frj,auij,stat_info)
#else
      SUBROUTINE rhs_force_ff_Newtonian_HuAdams_angular(this,&
           xi,xj,dij,vi,vj,numi,numj, pi,pj,&
           mi,mj,w,gradw,fi,fj,auij,stat_info)
#endif
        !----------------------------------------------------
        ! Subroutine : rhs_force_ff_Newtonian_HuAdams_angular
        !----------------------------------------------------
        !
        ! Purpose    : Implement the pressure, viscous and 
        !              thermal force with angular momentum 
        !              conservation.
        !
        ! Reference  : Hu and Adams
        !              Physics of Fluids 
        !              101702. 2006.
        !
        ! Remark     : Works only for incompressible flow.
        !
        ! Revisions  : V0.1 16.08.2010
        !
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
        ! 
        ! numi : number density of i
        ! numj : number density of j
        !----------------------------------------------------
        
        TYPE(Rhs), INTENT(INOUT)                :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xi
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xj
        REAL(MK), INTENT(IN)                    :: dij
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vi
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vj        
        REAL(MK), INTENT(IN)                    :: numi
        REAL(MK), INTENT(IN)                    :: numj
        REAL(MK), INTENT(IN)                    :: pi
        REAL(MK), INTENT(IN)                    :: pj
        REAL(MK), INTENT(IN)                    :: mi
        REAL(MK), INTENT(IN)                    :: mj
        REAL(MK), INTENT(IN)                    :: w
        REAL(MK), INTENT(IN)                    :: gradw
        REAL(MK), DIMENSION(:), INTENT(INOUT)   :: fi
        REAL(MK), DIMENSION(:), INTENT(INOUT)   :: fj
#ifdef __PARTICLES_FORCE_SEPARATE
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: fpi
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: fpj
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: fvi
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: fvj
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: fri
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: frj     
#endif
        REAL(MK), INTENT(OUT), OPTIONAL         :: auij
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_dim
        REAL(MK)                        :: dt
        REAL(MK)                        :: eta
        REAL(MK)                        :: eta_coef
        REAL(MK)                        :: kt
        REAL(MK), DIMENSION(3)          :: eij
        REAL(MK), DIMENSION(3)          :: vij
        REAL(MK)                        :: ev
        REAL(MK)                        :: f_c
        REAL(MK)                        :: f_d
        REAL(MK)                        :: f_r
        REAL(MK)                        :: pij
        REAL(MK)                        :: Aij
        REAL(MK)                        :: eps
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        fi(:) = 0.0_MK
        fj(:) = 0.0_MK
        
#ifdef __PARTICLES_FORCE_SEPARATE
        IF ( PRESENT(fpi) .AND. PRESENT (fpj) .AND. &
             PRESENT(fvi) .AND. PRESENT (fvj) .AND. &
             PRESENT(fri) .AND. PRESENT (frj) ) THEN
           
           fpi(:) = 0.0_MK
           fpj(:) = 0.0_MK
           fvi(:) = 0.0_MK
           fvj(:) = 0.0_MK
           fri(:) = 0.0_MK
           frj(:) = 0.0_MK 

        END IF
#endif
        num_dim = this%num_dim
        
        dt  = this%dt
        eta = this%eta
        eta_coef = this%eta_coef
        kt  = this%kt
        
        !----------------------------------------------------
        ! Check if any denominator is non-positive.
        !----------------------------------------------------
        
        IF ( mi < mcf_machine_zero .OR. &
             mj < mcf_machine_zero .OR. &
             numi < mcf_machine_zero .OR. &
             numj < mcf_machine_zero ) THEN

           PRINT *,"xi,xj,dij : ",  xi,xj,dij
           PRINT *,"vi,vj : ", vi,vj
           PRINT *,"numi,numj : ", numi,numj
           PRINT *,"pi,pj : ", pi,pj
           PRINT *,"mi,mj : ", mi,mj
           PRINT *,"w, gradw : ", w, gradw          
           
           stat_info = -1
           GOTO 9999
           
        END IF
        
        !----------------------------------------------------
        ! Calculate normalized vector pointing from j to i.
        !----------------------------------------------------
        
        eij(1:num_dim) = (xi(1:num_dim)  - xj(1:num_dim)) / dij
        
        !----------------------------------------------------
        ! Calculate the velocity vector pointing from j to i.
        !----------------------------------------------------
        
        vij(1:num_dim)  = vi(1:num_dim) - vj(1:num_dim)
        
        !----------------------------------------------------
        ! Dot product of eij and vij.
        !----------------------------------------------------
        
        ev = DOT_PRODUCT(eij(1:num_dim), vij(1:num_dim))

        pij = ( pi + pj ) / 2.0_MK
        Aij = ( 1.0_MK/numi**2.0_MK+ 1.0_MK/numj**2.0_MK ) * gradw
        
        !----------------------------------------------------
	! Calculate conservative force
	! per unit mass (from pressure).
  	!----------------------------------------------------
        
        f_c = Aij * pij
        
        fi(1:num_dim) = fi(1:num_dim) - &
             f_c * eij(1:num_dim) / mi
        fj(1:num_dim) = fj(1:num_dim) + &
             f_c * eij(1:num_dim) / mj
        
        !----------------------------------------------------
        ! Calculate potential energy, if needed.
        !----------------------------------------------------
        
        IF(PRESENT(auij))THEN
           auij = 0.5_MK * f_c * ev / mi
        END IF
        
        !----------------------------------------------------
        ! Calculate dissipative force
	! per unit mass(from viscosity).
  	!----------------------------------------------------
        
        f_d = eta_coef* eta * Aij * ev / dij
        
        fi(1:num_dim) = fi(1:num_dim) + &
             f_d  * eij(1:num_dim) / mi
        fj(1:num_dim) = fj(1:num_dim) - &
             f_d  * eij(1:num_dim) / mj
        
#ifdef __PARTICLES_FORCE_SEPARATE
        IF ( PRESENT(fpi) .AND. PRESENT (fpj) .AND. &
             PRESENT(fvi) .AND. PRESENT (fvj) ) THEN
             
           fpi(1:num_dim) = -f_c*eij(1:num_dim)/ mi
           fpj(1:num_dim) =  f_c*eij(1:num_dim)/ mj
           fvi(1:num_dim) =  f_d*eij(1:num_dim)/ mi
           fvj(1:num_dim) = -f_d*eij(1:num_dim)/ mj
           
        END IF
#endif
        !----------------------------------------------------
        ! Calculate random force
	! per unit mass(from thermal noise),
        ! if kt is above zero.
  	!----------------------------------------------------
        
        IF ( this%Brownian .AND. &
             kt >  mcf_machine_zero ) THEN
           
           eps = random_random(this%random,stat_info_sub)
           
           f_r = eps* SQRT(2.0_MK * kt * eta) * &
                SQRT(-eta_coef*Aij/dij) / SQRT(dt)
           
           fi(1:num_dim)  = fi(1:num_dim) + &
                f_r * eij(1:num_dim) / mi 
           fj(1:num_dim)  = fj(1:num_dim) - &
                f_r * eij(1:num_dim) / mj 
           
#ifdef __PARTICLES_FORCE_SEPARATE           
           IF ( PRESENT(fri) .AND. PRESENT (frj) ) THEN
              
              fri(1:num_dim) =  f_r*eij(1:num_dim)/ mi
              frj(1:num_dim) = -f_r*eij(1:num_dim)/ mj

           END IF
#endif
  
        END IF ! Brownian and kt > 0
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE rhs_force_ff_Newtonian_HuAdams_angular
      
      
