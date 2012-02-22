      SUBROUTINE statistic_compute_statistic(this,d_particles,&
           extra_k,extra_mom,stat_info)
        !----------------------------------------------------
        ! Subroutine  : statistic_compute_statistic
        !----------------------------------------------------
        ! Purpose     : Computes statistics, such as :
        !               total kinetic energy, total momentum.
        !               If there is colloid, add in the
        !               extra kinetic energy, extra momentum
        !               from colloids.
        !
        ! Routines    :
        !
        ! Remarks     :
        !
        ! References  :
        !
        ! Revisions   : V0.1 15.03 2009, original version.
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
        !
        ! this          :  Statistic object.
        ! d_particles   :  Particles object.
        ! extra_k       :  additional kinetic energy,
        !                  can be from colloidal particles.
        ! extra_mom     :  additional momentum,
        !                  can be from colloidal particles.
        ! stat_info     :  status information.
        !----------------------------------------------------
        
        TYPE(Statistic), INTENT(INOUT)          :: this
        TYPE(Particles), INTENT(IN)             :: d_particles
        REAL(MK),INTENT(IN), OPTIONAL           :: extra_k
        REAL(MK),DIMENSION(:),INTENT(IN),OPTIONAL :: extra_mom
        INTEGER, INTENT(OUT)                    :: stat_info

        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        TYPE(Control), POINTER                  :: ctrl
        LOGICAL                                 :: Brownian
        LOGICAL                                 :: stress_tensor
        LOGICAL                                 :: p_energy
        INTEGER                                 :: dim, dim2

        !----------------------------------------------------
        ! Particles variables,  
        ! position, velocity, mass and species ID.
        !----------------------------------------------------

        INTEGER                                 :: num_part
        REAL(MK), DIMENSION(:,:), POINTER       :: x
        REAL(MK), DIMENSION(:,:), POINTER       :: v
        REAL(MK), DIMENSION(:), POINTER         :: m
        INTEGER, DIMENSION(:), POINTER          :: sid
        REAL(MK), DIMENSION(:,:), POINTER       :: s
        REAL(MK), DIMENSION(:,:), POINTER       :: sp
        REAL(MK), DIMENSION(:,:), POINTER       :: sv
        REAL(MK), DIMENSION(:,:), POINTER       :: sr
        REAL(MK), DIMENSION(:), POINTER         :: u
        
        !----------------------------------------------------
        ! A technique object pointer.
        !----------------------------------------------------
        
        TYPE(Technique), POINTER                :: tech
        INTEGER                                 :: comm
        INTEGER                                 :: MPI_PREC
        
        !----------------------------------------------------
        ! Total values of kinetic energy, momentum.        
        !----------------------------------------------------
        
        REAL(MK)                                :: k_energy_tot
        REAL(MK), DIMENSION(3)                  :: momentum_tot
        REAL(MK)                                :: v2
        REAL(MK), DIMENSION(:), ALLOCATABLE     :: stress_tot
        REAL(MK), DIMENSION(:), ALLOCATABLE     :: stress_tot_p
        REAL(MK), DIMENSION(:), ALLOCATABLE     :: stress_tot_v
        REAL(MK), DIMENSION(:), ALLOCATABLE     :: stress_tot_r
        REAL(MK)                                :: p_energy_tot
        
        !----------------------------------------------------
        ! Integer counters and indices.
        !----------------------------------------------------
        
        INTEGER                                 :: i,j
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(ctrl)
        dim  = this%num_dim
        dim2 = dim**2
        
        NULLIFY(x)
        NULLIFY(v)
        NULLIFY(m)
        NULLIFY(sid)
        
        NULLIFY(s)
#ifdef __PARTICLES_STRESS_SEPARATE
        NULLIFY(sv)
        NULLIFY(sp)
        NULLIFY(sr)
#endif
        
        NULLIFY(u)
        
        NULLIFY(tech)
        
        this%k_energy    = 0.0_MK
        this%momentum(:) = 0.0_MK
        this%p_energy    = 0.0_MK
        
        k_energy_tot    = 0.0_MK
        momentum_tot(:) = 0.0_MK
        p_energy_tot    = 0.0_MK

        !----------------------------------------------------
        ! Check if dimension of additional momentum match.
        !----------------------------------------------------
        
        IF( PRESENT(extra_mom) .AND. &
             SIZE(extra_mom,1) /= dim) THEN
           PRINT *, "statistic_compute_statistic : ", &
                "Extra momentum has different dimension !"
           stat_info = -1
           GOTO 9999           
        END IF
        
        !----------------------------------------------------
        ! Get control parameters.
        !----------------------------------------------------
   
        CALL particles_get_ctrl(d_particles,ctrl,stat_info_sub)
        Brownian = &
             control_get_Brownian(ctrl,stat_info_sub)     
        stress_tensor = &
             control_get_stress_tensor(ctrl,stat_info_sub)
        p_energy = &
             control_get_p_energy(ctrl,stat_info_sub)
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        IF ( stress_tensor ) THEN
           
           this%stress(:)  = 0.0_MK
           ALLOCATE(stress_tot(dim2))
           stress_tot(:)   = 0.0_MK 
           
#ifdef __PARTICLES_STRESS_SEPARATE
           this%stress_p(:)  = 0.0_MK
           ALLOCATE(stress_tot_p(dim2))
           stress_tot_p(:)   = 0.0_MK 
           this%stress_v(:)  = 0.0_MK
           ALLOCATE(stress_tot_v(dim2))
           stress_tot_v(:)   = 0.0_MK 
           
           IF ( Brownian ) THEN
              this%stress_r(:)  = 0.0_MK
              ALLOCATE(stress_tot_r(dim2))
              stress_tot_r(:)   = 0.0_MK         
           END IF
#endif
           

        END IF
        
        !----------------------------------------------------
        ! Get rank of this process.
        !----------------------------------------------------
        
        CALL particles_get_tech(d_particles,tech,stat_info_sub)
        
        !----------------------------------------------------
        ! Get number of real particles on this process.
        !----------------------------------------------------
        
        num_part = &
             particles_get_num_part_real(d_particles,stat_info_sub)
        
        !----------------------------------------------------
        ! Get position, velocity, mass, speciesID of particles.
        !----------------------------------------------------
        
        CALL particles_get_x(d_particles,x,num_part,stat_info_sub)
        CALL particles_get_v(d_particles,v,num_part,stat_info_sub)
        CALL particles_get_m(d_particles,m,num_part,stat_info_sub)
        CALL particles_get_sid(d_particles,sid,num_part,stat_info_sub)

        IF ( stress_tensor) THEN
           
           CALL particles_get_s(d_particles,s,num_part,stat_info_sub)
#ifdef __PARTICLES_STRESS_SEPARATE
           CALL particles_get_sp(d_particles,sp,num_part,stat_info_sub)
           CALL particles_get_sv(d_particles,sv,num_part,stat_info_sub)
           IF ( Brownian ) THEN
              CALL particles_get_sr(d_particles,sr,num_part,stat_info_sub)
           END IF
#endif
           
        END IF
     
        !----------------------------------------------------
        ! Kinetic engergy : K = 0.5 * m * v**2; 
        ! Momentum        : M = m * v; 
        !----------------------------------------------------
        
        DO j = 1, num_part
           
           !-------------------------------------------------
           ! Count only fluid particles here 
           ! when calculate total kinetic energy and momentum.
           !-------------------------------------------------
           
           IF( sid(j) == 0 ) THEN
              
              v2 = 0.0_MK
              
              DO i = 1, dim
                 
                 v2 = v2 + v(i,j)**2
                 
              END DO
              
              this%k_energy  = this%k_energy + &
                   0.5_MK * m(j) * v2
              
              this%momentum(1:dim) = &
                   this%momentum(1:dim) + &
                   m(j) * v(1:dim,j)
              
           END IF ! sid(j) == 0
           
           !-------------------------------------------------
           ! Count solvent particles for total stress tensor
           !-------------------------------------------------
           
           IF ( stress_tensor ) THEN
              
              IF ( sid(j) == 0 ) THEN
                 
                 this%stress(1:dim2) = &
                      this%stress(1:dim2) + s(1:dim2,j)
                 
#ifdef __PARTICLES_STRESS_SEPARATE
                 this%stress_p(1:dim2) = &
                      this%stress_p(1:dim2) + sp(1:dim2,j)
                 this%stress_v(1:dim2) = &
                      this%stress_v(1:dim2) + sv(1:dim2,j)
                 
                 IF ( Brownian ) THEN
                    this%stress_r(1:dim2) = &
                         this%stress_r(1:dim2) + sr(1:dim2,j)
                 END IF
#endif
                 
              END IF ! sid
              
           END IF ! stress_tensor
           
        END DO ! j = 1, num_part
        
        
        IF( p_energy ) THEN
           
           CALL particles_get_u(d_particles,u,num_part,stat_info_sub)
           
           DO j = 1, num_part
              
              !----------------------------------------------
              ! Count only fluid particles.
              !----------------------------------------------
              
              IF( sid(j) == 0 ) THEN
                 
                 this%p_energy = &
                      this%p_energy + m(j)*u(j)
                 
              END IF
              
           END DO
           
        END IF

#if 0        
        !----------------------------------------------------
        ! Commented out this area, as stress on each particle
        ! has been divided by 2 already.
        !----------------------------------------------------
        !----------------------------------------------------
        ! As each pair particle is counted twice, 
        ! we divide 2 here.
        ! It is done before MPI summation to have consistency
        ! for MPI and non-MPI runs.
        !----------------------------------------------------
        
        IF ( stress_tensor ) THEN
           
           this%stress(1:dim2) = this%stress(1:dim2) / 2.0_MK
           
#ifdef __PARTICLES_STRESS_SEPARATE
           this%stress_p(1:dim2) = this%stress_p(1:dim2) / 2.0_MK
           this%stress_v(1:dim2) = this%stress_v(1:dim2) / 2.0_MK
           IF ( Brownian ) THEN
              this%stress_r(1:dim2) = this%stress_r(1:dim2) / 2.0_MK
           END IF
#endif
           
        END IF
#endif  
        
        !----------------------------------------------------
        ! Calculation in the context of MPI.
        !----------------------------------------------------
        
#ifdef __MPI
        
        comm     = technique_get_comm(tech,stat_info_sub)
        MPI_PREC = technique_get_MPI_PREC(tech,stat_info_sub)
        
        !----------------------------------------------------
        ! Sum up the kinetic energy from all processes,
        ! broadcast the result.
        !----------------------------------------------------
        
        CALL MPI_ALLREDUCE (this%k_energy,k_energy_tot, &
             1,MPI_PREC,MPI_SUM,comm,stat_info_sub) 
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "statistic_compute_statistic : ", &
                "MPI_ALLREDUCE() for kinetic energy has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Sum up the momentum from all the processes,
        ! broadcast the result.
        !----------------------------------------------------
        
        CALL MPI_ALLREDUCE (this%momentum(1:dim),  &
             momentum_tot(1:dim),dim,MPI_PREC, &
             MPI_SUM,comm,stat_info_sub) 
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "statistic_compute_statistic : ", &
                "MPI_ALLREDUCE() for momentum has problem !"
           stat_info = -1
           GOTO 9999
        END IF

        
        !----------------------------------------------------
        ! Sum up the stress tensor from all the processes,
        ! broadcast the result.
        !----------------------------------------------------
        
        IF ( stress_tensor ) THEN
           
           CALL MPI_ALLREDUCE (this%stress(1:dim2),  &
                stress_tot(1:dim2),dim2,MPI_PREC, &
                MPI_SUM,comm,stat_info_sub) 

#ifdef __PARTICLES_STRESS_SEPARATE
           CALL MPI_ALLREDUCE (this%stress_p(1:dim2),  &
                stress_tot_p(1:dim2),dim2,MPI_PREC, &
                MPI_SUM,comm,stat_info_sub) 
           CALL MPI_ALLREDUCE (this%stress_v(1:dim2),  &
                stress_tot_v(1:dim2),dim2,MPI_PREC, &
                MPI_SUM,comm,stat_info_sub) 
           IF ( Brownian ) THEN
              CALL MPI_ALLREDUCE (this%stress_r(1:dim2),  &
                   stress_tot_r(1:dim2),dim2,MPI_PREC, &
                   MPI_SUM,comm,stat_info_sub) 
           END IF
           
#endif
           
           IF( stat_info_sub /= 0 ) THEN
              PRINT *, "statistic_compute_statistic : ", &
                   "MPI_ALLREDUCE() for stress has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! All processes get total kinetic energy, 
        ! momentum and total flow rate.
        !----------------------------------------------------
        
        this%k_energy = k_energy_tot
        this%momentum(1:dim) = momentum_tot(1:dim)
        
        
        !----------------------------------------------------
        ! All processes get total stress tensor.
        !----------------------------------------------------
        
        IF ( stress_tensor ) THEN
           
           this%stress(1:dim2)  = stress_tot(1:dim2)
           
#ifdef __PARTICLES_STRESS_SEPARATE
           this%stress_p(1:dim2)  = stress_tot_p(1:dim2)
           this%stress_v(1:dim2)  = stress_tot_v(1:dim2)
           IF ( Brownian ) THEN
              this%stress_r(1:dim2)  = stress_tot_r(1:dim2)
           END IF
#endif
           
        END IF
        
        !----------------------------------------------------
        ! Sum up the potential energy from all processes,
        ! then broadcast the result.
        ! All processes get total potential energy.
        !----------------------------------------------------
        
        IF ( p_energy ) THEN
           
           CALL MPI_ALLREDUCE (this%p_energy,p_energy_tot, &
                1,MPI_PREC, MPI_SUM,comm,stat_info_sub) 
           
           IF( stat_info_sub /= 0 ) THEN
              PRINT *, "statistic_compute_statistic : ", &
                   "MPI_ALLREDUCE() for potential energy has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           this%p_energy = p_energy_tot
           
        END IF
        
#endif !__MPI
        
        !----------------------------------------------------
        ! Add extra kinetic energy and momentum.
        !----------------------------------------------------
        
        IF ( PRESENT(extra_k) ) THEN
           
           this%k_energy = this%k_energy + extra_k
           
        END IF
        
        IF ( PRESENT(extra_mom)) THEN
           
           this%momentum(1:dim) =  &
                this%momentum(1:dim) + extra_mom(1:dim)
           
        END IF
        
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Release the memory pointed by potiners.
        !----------------------------------------------------
        
        IF(ASSOCIATED(x)) THEN
           DEALLOCATE(x)
        END IF
        
        IF(ASSOCIATED(v)) THEN
           DEALLOCATE(v)
        END IF
        
        IF(ASSOCIATED(m)) THEN
           DEALLOCATE(m)
        END IF
        
        IF(ASSOCIATED(sid)) THEN
           DEALLOCATE(sid)
        END IF
        
        IF(ASSOCIATED(s)) THEN
           DEALLOCATE(s)
        END IF
        
#ifdef __PARTICLES_STRESS_SEPARATE
        IF(ASSOCIATED(sp)) THEN
           DEALLOCATE(sp)
        END IF
        IF(ASSOCIATED(sv)) THEN
           DEALLOCATE(sv)
        END IF
        IF(ASSOCIATED(sr)) THEN
           DEALLOCATE(sr)
        END IF
#endif
        
        IF(ASSOCIATED(u)) THEN
           DEALLOCATE(u)
        END IF
        
        
        RETURN          
        
      END SUBROUTINE statistic_compute_statistic
      
      
