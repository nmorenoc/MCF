      SUBROUTINE  particles_init_partial_inter(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_init_partial_inter
        !----------------------------------------------------
        !
        ! Purpose     : Create the rest quantities,
        !               besides the ones generated globally 
        !               already on root process.
        !                  
        ! Routines    :
        !
        ! Remarks     : Create the rest quntities,
        !               besides position, velocity,
        !               p_ID, s_ID for fluid particles 
        !               in sub-domains of each process.
        !
        !               Such as mass.
        !
        !               In case of non-Newtonian viscoelastic
        !               Oldroyd-B model fluid,
        !               conformation tensor has to be
        !               allocated.
        !
        !               for eigen-dynamics :
        !               eigenvalue and eigenvector have to
        !               be allocated.
        !               for evolution of conformation tensor :
        !               acceleration of conformation tensor
        !               will be allocated during calculation
        !               later.
        !
        !               We allocate potential energy array
        !               if needed.
        !
        ! References  :
        !
        ! Revisions   : V0.2 04.12 2009, check work flow
        !               and supply with more comments.
        !
        !               V0.1 30.07 2009, original version.
        !                 
        !----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        
        !----------------------------------------------------
        ! Arguments
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER,INTENT(OUT)	                :: stat_info
    	
	!----------------------------------------------------
    	! Local variables :
        !
        ! Newtonian : indicate if it is Newtonian fluid.
        ! p_energy  : indicate if potential energy needed.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        
        LOGICAL                                 :: Newtonian
        LOGICAL                                 :: p_energy
        INTEGER                                 :: num_dim,j
        LOGICAL                                 :: eigen_dynamics
        REAL(MK), DIMENSION(:), POINTER         :: eval
        REAL(MK), DIMENSION(:,:), POINTER       :: evec
        INTEGER                                 :: num_part_real
        INTEGER                                 :: i
        
        
        !----------------------------------------------------
    	! Initialization of variables.
    	!----------------------------------------------------
        
        stat_info = 0
        stat_info_sub = 0
        
        NULLIFY(eval)
        NULLIFY(evec)
        
        !----------------------------------------------------
        ! control variables.
        !----------------------------------------------------
        
        Newtonian = &
             control_get_Newtonian(this%ctrl,stat_info_sub)        
        p_energy  = &
             control_get_p_energy(this%ctrl,stat_info_sub)
        
        !----------------------------------------------------
        ! physics variables.
        !----------------------------------------------------
        
        num_dim       = &
             physics_get_num_dim(this%phys,stat_info_sub)
        
        num_part_real = this%num_part_real
        
        
        !----------------------------------------------------
        ! Allocate memory for mass.
        !----------------------------------------------------
        
        IF (ASSOCIATED(this%m)) THEN
           DEALLOCATE(this%m)
        END IF
        
        ALLOCATE(this%m(num_part_real),STAT=stat_info_sub)
        this%m(1:num_part_real) = 0.0_MK
        
        !----------------------------------------------------
        ! Allocate memory for conformation tensor
        ! for real particles, in case we are dealing
        ! with non-Newtonian fluid.
        !----------------------------------------------------
        
        IF( .NOT. Newtonian ) THEN
           
           IF (ASSOCIATED(this%ct))THEN
              DEALLOCATE(this%ct)
           END IF
           
           ALLOCATE(this%ct(num_dim**2,num_part_real), &
                STAT=stat_info_sub)
           
           IF(stat_info_sub /= 0) THEN
              PRINT *, "particles_init_partial_inter : ", &
                   "Allocating ct has problem !"
              stat_info = -1 
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Initialize ct as unit tensor.
           !-------------------------------------------------
           
           this%ct(:,:) = 0.0_MK
           
           DO i =1, num_dim
              this%ct(i+num_dim*(i-1),1:num_part_real) = 1.0_MK 
           END DO
           
           
           eigen_dynamics = &
                physics_get_eigen_dynamics(this%phys,stat_info_sub)
           
           !-------------------------------------------------
           ! For egenvector dynamics, we need egenvalues
           ! and egenvectors also.
           !-------------------------------------------------
           
           IF ( eigen_dynamics ) THEN
              
              IF (ASSOCIATED(this%eval))THEN
                 DEALLOCATE(this%eval)
              END IF
              
              IF (ASSOCIATED(this%evec))THEN
                 DEALLOCATE(this%evec)
              END IF
              
              ALLOCATE(this%eval(num_dim,num_part_real),&
                   STAT=stat_info_sub)
              
              ALLOCATE(this%evec(num_dim**2,num_part_real),&
                   STAT=stat_info_sub)
              
              CALL physics_get_eval(this%phys,eval,stat_info_sub)
              CALL physics_get_evec(this%phys,evec,stat_info_sub)
              
              !----------------------------------------------
              ! Set the initial values for eigenvalues and
              ! eigenvectors
              !----------------------------------------------
              
              DO i = 1, num_dim
                 this%eval(i,1:num_part_real) = eval(i)
              END DO
              
              DO j = 1, num_dim
                 DO i = 1, num_dim
                    this%evec(i+num_dim*(j-1),1:num_part_real) = evec(i,j)
                 END DO
              END DO
              
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! u  : potential energy.
        !      It is needed for testing code only,
        !----------------------------------------------------
        
        IF( p_energy ) THEN
           
           IF (ASSOCIATED(this%u)) THEN
              DEALLOCATE(this%u,STAT=stat_info_sub)
           END IF
           
           ALLOCATE(this%u(num_part_real), &
                STAT=stat_info_sub)
           
           IF(stat_info_sub /= 0) THEN
              PRINT *, &
                   "particles_init_partial_inter : ",&
                   "Allocating u has problem !"
              stat_info = -1 
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Initialize potential energy per 
           ! unit mass for each particle as 0.5.
           !-------------------------------------------------
           
           this%u(1:num_part_real) = 0.5_MK           
           
        END IF
        
        
9999	CONTINUE
        
        IF( ASSOCIATED(eval)) THEN
           DEALLOCATE(eval)
        END IF
        
        IF( ASSOCIATED(evec)) THEN
           DEALLOCATE(evec)
        END IF
        
	RETURN
 
      END SUBROUTINE particles_init_partial_inter

      
     
