!----------------------------------------------------------------
! This contains the routines which calculate the force and 
! return two accelerations between two colloidal boundary
! particle Newtonian particles, 
! accelerations only from conservative(pressure) force 
!----------------------------------------------------------------

      SUBROUTINE rhs_force_cc_Newtonian(this,&
           xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&      
           mi,mj,w,gradw,fi,fj,auij,stat_info)
        !----------------------------------------------------
        ! Subroutine  :  Return pair-wise force per unit
        !                mass, i.e., acceleration between
        !                two Newtonian fluid particles.
        !                Currently, even we have colloids
        !                boundary particles, 
        !                fluid-boundary particles force
        !                calculation is treated same way.
        ! 
        !  Revision   :               
        !               V0.1 05.05.2010, original version
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
        !  Arguments
        !  xi        : particle i 's position
        !  xj        : particle j 's position
        !  dij       : distance between i,j
        !  vi        : i's velocity
        !  vj        : j's velocity
        !  rhoi      : i's density
        !  rhoj      : j's density
        !  pi        : i's pressure
        !  pj        : j's pressure
        !  mi        : i's mass
        !  mj        : j's mass
        !  w         : kernel value
        !  gradw     : kernel gradient value
        !  fi        : force acting on i (per unit mass)
        !  fj        : force acting on j (per unit mass)
        !  auij      : potential energy between i,j
        !  stat_info : return flag of status
        !----------------------------------------------------
        
        TYPE(Rhs), INTENT(IN)                   :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xi
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xj
        REAL(MK), INTENT(IN)                    :: dij
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vi
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vj
        REAL(MK), INTENT(IN)                    :: rhoi
        REAL(MK), INTENT(IN)                    :: rhoj
        REAL(MK), INTENT(IN)                    :: pi
        REAL(MK), INTENT(IN)                    :: pj
        REAL(MK), INTENT(IN)                    :: mi
        REAL(MK), INTENT(IN)                    :: mj
        REAL(MK), INTENT(IN)                    :: w
        REAL(MK), INTENT(IN)                    :: gradw
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: fi
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: fj
        REAL(MK), INTENT(OUT), OPTIONAL         :: auij 
        INTEGER, INTENT(OUT)                    :: stat_info
        
        
        INTEGER                                 :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        fi(1:this%num_dim) = 0.0_MK
        fj(1:this%num_dim) = 0.0_MK
        
        IF(PRESENT(auij))THEN
           auij = 0.0_MK
        END IF
        
        SELECT CASE (this%rhs_force_type)
           
        CASE (0) 
           
           PRINT *, "rhs_force_cc_Newtonian : ", &
                "0 not existing formulation !"
           stat_info = -1
           GOTO 9999
           
        CASE (1)
           
           PRINT *, "rhs_force_cc_Newtonian : ", &
                "1 not existing formulation !"
           stat_info = -1
           GOTO 9999
           
        CASE (2)
           
           PRINT *, "rhs_force_cc_Newtonian : ", &
                "2 not existing formulation !"
           stat_info = -1
           GOTO 9999
           
        CASE (3)
           
           CALL  rhs_force_cc_Newtonian_HuAdams(this,&
                xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&
                mi,mj,w,gradw,fi,fj,auij,stat_info_sub)
           
        CASE (4:) 
           
           PRINT *, "rhs_force_cc_Newtonian : ", & 
                "4-> Not existing formulation !"
           stat_info = -1
           GOTO 9999
           
        END SELECT
        
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "rhs_force_cc_Newtonian :", &
                "Formulation has some problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE rhs_force_cc_Newtonian
      
      
      
      SUBROUTINE rhs_force_cc_Newtonian_HuAdams(this,&
           xi,xj,dij,vi,vj,numi,numj, pi,pj,&
           mi,mj,w,gradw,fi,fj,auij,stat_info)
        !----------------------------------------------------
        ! Subroutine : rhs_force_ff_Newtonian_HuAdams
        !----------------------------------------------------
        !
        ! Purpose    : Implement the pressure, viscous and 
        !              thermal force from Hu and Adams 2006 
        !              Journal of Computational Physics.
        !
        ! Reference  : X.Y.Hu and N.A.Adams Journal of 
        !              Computational Physics 213 (2006)
        !              844-861.
        !
        ! Remark     :
        !
        ! Revisions  : V0.1 05.05.2010
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
        
        TYPE(Rhs), INTENT(IN)                   :: this
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
        REAL(MK), INTENT(OUT), OPTIONAL         :: auij
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_dim
        REAL(MK)                        :: eta
        REAL(MK)                        :: kt
        REAL(MK), DIMENSION(3)          :: eij
        REAL(MK), DIMENSION(3)          :: vij
        REAL(MK)                        :: ev
        REAL(MK)                        :: f_c
        INTEGER                         :: i
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        num_dim = this%num_dim
        
        eta = physics_get_eta(this%phys,stat_info_sub)     
        kt  = physics_get_kt(this%phys,stat_info_sub)
        
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
        
        ev = 0.0_MK
        DO i = 1, num_dim
           ev = ev + eij(i)*vij(i)
        END DO
        
        !----------------------------------------------------
	! Calculate conservative force
	! per unit mass (from pressure).
  	!----------------------------------------------------
        
        f_c = ( pi/(numi**2.0_MK)+ pj/(numj**2.0_MK)) * gradw
        
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
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE rhs_force_cc_Newtonian_HuAdams
      
      
