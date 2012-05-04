      SUBROUTINE rhs_force_ff_non_Newtonian_HuAdams_angular(this,&
           xi,xj,dij,vi,vj,numi,numj,pi,pj,&
           mi,mj,w,gradw,fi,fj,stat_info)
        !----------------------------------------------------
        ! Subroutine : rhs_force_ff_non_Newtonian_HuAdams_angular
        !----------------------------------------------------
        !----------------------------------------------------
        !
        ! Purpose    : Implement the pressure, viscous and 
        !              thermal force from
        !              Hu and Adams Physics of fluids 2006,
        !              which conserve angular momentum.
        !              In addition, With extra pressure tensor,
        !              to handle a viscoelastic fluid.
        !
        ! Reference  : Hu&Adams, Physics of Fluids 101702. 2006.
        !              Vazquez, Eller, Espanol, Phys.Rev.E.2009. 
        !
        ! Remark     : Works only for incompressible flow and
        !              deterministic flow. i.e., the Brownian
        !              part is not implemented.
        !
        ! Revisions  : V0.1 04.02.2012
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
     
        TYPE(Rhs), INTENT(INOUT)                :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xi
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xj
        REAL(MK), INTENT(IN)                    :: dij
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vi
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vj
        !number density 
        REAL(MK), INTENT(IN)                    :: numi
        REAL(MK), INTENT(IN)                    :: numj
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: pi
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: pj
        REAL(MK), INTENT(IN)                    :: mi
        REAL(MK), INTENT(IN)                    :: mj
        REAL(MK), INTENT(IN)                    :: w
        REAL(MK), INTENT(IN)                    :: gradw
        REAL(MK), DIMENSION(:), INTENT(INOUT)   :: fi
        REAL(MK), DIMENSION(:), INTENT(INOUT)   :: fj
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
        REAL(MK)                        :: pij        
        REAL(MK), DIMENSION(3)          :: eij
        REAL(MK), DIMENSION(3)          :: vij
        REAL(MK)                        :: ev
        REAL(MK), DIMENSION(3)          :: f_c
        REAL(MK)                        :: f_d
        INTEGER                         :: i,j
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        fi(:) = 0.0_MK
        fj(:) = 0.0_MK
        
        num_dim = this%num_dim
        
        dt        = this%dt
        eta       = this%eta
        eta_coef  = this%eta_coef
        kt        = this%kt
        
        !----------------------------------------------------
        ! Check if any denominator is non-positive.
        !----------------------------------------------------
        
        IF ( dij < mcf_machine_zero .OR. &
             mi < mcf_machine_zero .OR. &
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
        
        eij(1:num_dim) = (xi(1:num_dim) - xj(1:num_dim)) / dij
        
        !----------------------------------------------------
        ! Calculate velocity vector pointing from j to i.
        !----------------------------------------------------
        
        vij(1:num_dim) = vi(1:num_dim) - vj(1:num_dim)
        
        !----------------------------------------------------
        ! Dot product of eij and vij.
        !----------------------------------------------------
        
        ev = DOT_PRODUCT(eij(1:num_dim), vij(1:num_dim))

        !----------------------------------------------------
	! Calculate the conservative force
	! per unit mass(from pressure).
  	!----------------------------------------------------
        
        f_c(1:num_dim) = 0.0_MK

  	!----------------------------------------------------
        ! According to Hu&Adams, Physics of Fluids 2006,
        ! the pressure between i,j is averaged before
        ! calculating pressure force.
        ! Therefore, we average the diagonal pressure tensor
        ! between i and j.
  	!----------------------------------------------------
        
        DO j = 1, num_dim
           DO i = 1, num_dim
              IF ( i == j ) THEN
                 pij = ( pi(i,j) + pj(i,j) ) / 2.0_MK
                 f_c(j) = f_c(j) + &
                      pij * eij(i) / (numi**2.0_MK) +  &
                      pij * eij(i) /( numj**2.0_MK)
              ELSE
                 f_c(j) = f_c(j) + &
                      pi(i,j) * eij(i) / (numi**2.0_MK) +  &
                      pj(i,j) * eij(i) /( numj**2.0_MK)
              END IF
           END DO
        END DO
        
        fi(1:num_dim) = fi(1:num_dim) - &
             f_c(1:num_dim) * gradw / mi
        fj(1:num_dim) = fj(1:num_dim) + &
             f_c(1:num_dim) * gradw / mj
        
        !----------------------------------------------------
        ! Calculate the dissipative force
	! per unit mass(from viscosity).
  	!----------------------------------------------------
        
        f_d =  eta_coef * eta * (1.0_MK/numi**2.0_MK + &
             1.0_MK/numj**2.0_MK) * gradw * ev / dij
        
        fi(1:num_dim) = fi(1:num_dim) + &
             f_d  * eij(1:num_dim) / mi
        fj(1:num_dim) = fj(1:num_dim) - &
             f_d  * eij(1:num_dim) / mj
        
        !----------------------------------------------------
        ! Calculate random force
	! per unit mass(from thermal noise),
        ! if kt is above zero.
  	!----------------------------------------------------
        
        IF ( this%Brownian .AND. &
             kt >  mcf_machine_zero ) THEN
           
           PRINT *, "rhs_force_ff_non_Newtonian_HuAdams_agnular: ",&
                "Brownian force is not available !"
           stat_info = -1
           GOTO 9999
           
        END IF ! Brownian and kt > 0
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE rhs_force_ff_non_Newtonian_HuAdams_angular
      
