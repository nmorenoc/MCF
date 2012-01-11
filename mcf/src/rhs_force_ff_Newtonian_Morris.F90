!----------------------------------------------------------------
! This contains the routines which calculate the force and 
! return two accelerations between two Newtonian particles, 
! accelerations from i.e. conservative(pressure) force 
! dissipative(viscous) force and random(thermal noise) force.
! using Morris et al 1997 formulation.
!----------------------------------------------------------------
     
      SUBROUTINE rhs_force_ff_Newtonian_Morris(this,&
           xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&
           mi,mj,w,gradw,fi,fj,auij,stat_info)
        !----------------------------------------------------
        ! Implementing pressure and viscous force from
        ! Morris et al. 1997 JCP paper.
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
        REAL(MK), DIMENSION(:), INTENT(INOUT)   :: fi
        REAL(MK), DIMENSION(:), INTENT(INOUT)   :: fj
        REAL(MK), INTENT(OUT),OPTIONAL          :: auij
        
        INTEGER, INTENT(OUT)                    :: stat_info
        
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_dim,i
        REAL(MK), DIMENSION(3)          :: eij
        REAL(MK)                        :: ev
        REAL(MK), DIMENSION(3)          :: vij
        REAL(MK)                        :: f_c,f_d
        REAL(MK)                        :: eta
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        fi(:) = 0.0_MK
        fj(:) = 0.0_MK
        
        num_dim = this%num_dim
        eta     = physics_get_eta(this%phys,stat_info_sub)                
        
        
        !----------------------------------------------------
        ! Check if any denominator is non-positive.
        !----------------------------------------------------
        
        IF ( dij < mcf_machine_zero .OR. &
             rhoi < mcf_machine_zero .OR. &
             rhoj < mcf_machine_zero ) THEN 
           
           PRINT *,"xi,xj,dij : ",  xi,xj,dij
           PRINT *,"vi,vj : ", vi,vj
           PRINT *,"rhoi,rhoj : ", rhoi,rhoj
           PRINT *,"pi,pj : ", pi,pj
           PRINT *,"mi,mj : ", mi,mj
           PRINT *,"w, gradw : ", w, gradw
           
           PRINT *, "rhs_force_ff_Newtonian_Morris : ",& 
                "Divided by non-positive !" 
           stat_info = -1
           GOTO 9999           
        END IF
        
        !----------------------------------------------------
	! Calculate normalized vector pointing from j to i
 	!----------------------------------------------------
        
        eij(1:num_dim) = (xi(1:num_dim) - xj(1:num_dim)) / dij
        
 	!----------------------------------------------------
	! Calculate velocity vector pointing from j to i
 	!----------------------------------------------------
        
        vij(1:num_dim)  = vi(1:num_dim) - vj(1:num_dim)
        
        !----------------------------------------------------
        ! Dot product of eij and vij.
 	!----------------------------------------------------
        
        ev = 0.0_MK        
        DO i = 1,num_dim
           ev = ev + eij(i) * vij(i)
        END DO
        
  	!------------------------------------------
	! Calculate the conservative force
	! per unit mass (from pressure).
        ! Accleration of thermal energy.
  	!------------------------------------------
        
        f_c = ( pi/(rhoi**2.0_MK) + pj/(rhoj**2.0_MK)) * gradw
        
        fi(1:num_dim) = fi(1:num_dim) - &
             mj * f_c * eij(1:num_dim) 
        fj(1:num_dim) = fj(1:num_dim) + &
             mi * f_c * eij(1:num_dim) 
        
        IF ( PRESENT(auij) ) THEN
           auij  = 0.5_MK * mj * f_c * ev 
        END IF
        
  	!------------------------------------------
	! Calculate the dissipative force
	! per unit mass(from viscosity)
  	!------------------------------------------
        
        f_d = 2.0_MK * eta * gradw / rhoi / rhoj / dij
        
        fi(1:num_dim)  = fi(1:num_dim) + &
             mj * f_d  * vij(1:num_dim)  
        fj(1:num_dim)  = fj(1:num_dim) - &
             mi * f_d  * vij(1:num_dim)  
        
9999 	CONTINUE
        
        RETURN
        
      END SUBROUTINE rhs_force_ff_Newtonian_Morris

