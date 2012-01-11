      SUBROUTINE particles_integrate_evec(this,&
           num,dt,lambda,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_integrate_evec
        !----------------------------------------------------
        !
        ! Purpose     : Integrate egenvectors of particles
        !               with required accuracy, in case of
        !               egenvector dynamics.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revision    :  V0.1 04.08.2009, original version.
        !
        !----------------------------------------------------
        !  Author    : Xin Bian
        !  Contact   : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        !  Arguments
        !
        !  this       : an object of Particles Class.
        !  num        : number of particles needed to be updated,
        !               i.e. first num particles in this%x 
        !               are operated.
        !  dt         : time step.
        !  lambda     : coefficient required.
        !  stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        REAL(MK), INTENT(IN)                    :: dt
        REAL(MK), INTENT(IN)                    :: lambda
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim, dim2,i, j, k
        REAL(MK)                                :: len, evec_tolerance
        
        !----------------------------------------------------
        !  Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        IF( num > this%num_part_real) THEN
           PRINT *, "particles_integrate_evec : ", &
                "num > num_part_real, wrong !"
           stat_info = -1
           GOTO 9999      
        END IF
        
        dim = physics_get_num_dim(this%phys,stat_info_sub)
        evec_tolerance = &
             physics_get_evec_tolerance(this%phys,stat_info_sub)
        
        dim2 = dim**2
        
        !----------------------------------------------------
        ! Integrate eigenvectors.
        !----------------------------------------------------
        
        DO i = 1, num
           
           this%evec(1:dim2,i) = &
                this%evec(1:dim2,i) + &
                this%aevec(1:dim2,i) * dt * lambda
           
           DO j = 1, dim

              k = (j-1) * dim + 1
              
              len = SQRT(DOT_PRODUCT(this%evec(k:j*dim,i),this%evec(k:j*dim,i)))
              
              IF (ABS(1.0_MK - len) > evec_tolerance) THEN
                 
                 this%evec(k:j*dim,i) = this%evec(k:j*dim,i) / len
                 
              END IF
              
           END DO ! j
           
        END DO ! i 
        
9999    CONTINUE      
        
        RETURN
        
      END SUBROUTINE particles_integrate_evec
      
      
