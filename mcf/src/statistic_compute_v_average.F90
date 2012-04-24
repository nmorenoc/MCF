      SUBROUTINE statistic_compute_v_average(this,d_particles,stat_info)
        !----------------------------------------------------
        ! Subroutine  : statistic_compute_v_average
        !----------------------------------------------------
        ! Purpose     : Computes average flow velocity
        !               in certain width at lower side of
        !               the box.
        !             
        ! Routines    :
        !
        ! Remarks     :
        !
        ! References  :
        !
        ! Revisions   : V0.1 09.12 2009, original version.
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
        ! stat_info     :  status information.
        !----------------------------------------------------
        
        TYPE(Statistic), INTENT(INOUT)          :: this
        TYPE(Particles), INTENT(IN)             :: d_particles
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        
        !----------------------------------------------------
        ! phys     : a physics object pointer.
        ! num_dim  : number of dimension.
        ! fd       : flow direction.
        ! fw         The width rate of domain to calculate 
        !            flow rate.
        !----------------------------------------------------
        
        TYPE(Physics), POINTER                  :: phys
        INTEGER                                 :: num_dim, j
        REAL(MK), DIMENSION(:), POINTER         :: max_phys
        REAL(MK), DIMENSION(:), POINTER         :: min_phys
        INTEGER                                 :: fd
        REAl(MK)                                :: fw
        
        !----------------------------------------------------
        ! Particles variables,  
        ! position, velocity, mass and species ID.
        !----------------------------------------------------
        
        INTEGER                                 :: num_part
        REAL(MK), DIMENSION(:,:), POINTER       :: x
        REAL(MK), DIMENSION(:,:), POINTER       :: v
        REAL(MK), DIMENSION(:), POINTER         :: m
        INTEGER, DIMENSION(:), POINTER          :: sid
        
        !----------------------------------------------------
        ! A technique object pointer.
        !----------------------------------------------------
         
        TYPE(Technique), POINTER                :: tech
        INTEGER                                 :: comm
        INTEGER                                 :: MPI_PREC
        
        
        !----------------------------------------------------
        ! Auxiliary parameters :
        !
        ! Integer counters and indices.
        !----------------------------------------------------
        
        REAL(MK), DIMENSION(3)                  :: v_aver_tot
        REAL(MK)                                :: v_mass
        REAL(MK)                                :: v_mass_tot        
        

        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(phys)
        num_dim = this%num_dim
        NULLIFY(max_phys)
        NULLIFY(min_phys)      
        
        
        NULLIFY(x)
        NULLIFY(v)
        NULLIFY(m)
        NULLIFY(sid)
        
        NULLIFY(tech)
        
        v_aver_tot(:)   = 0.0_MK
        this%v_aver(:)  = 0.0_MK
        v_mass          = 0.0_MK
        v_mass_tot      = 0.0_MK
        
        !----------------------------------------------------
        ! Get physics parameters.
        !----------------------------------------------------
        
        CALL particles_get_phys(d_particles,phys,stat_info_sub)
        CALL physics_get_max_phys(phys,max_phys,stat_info_sub)
        CALL physics_get_min_phys(phys,min_phys,stat_info_sub)
        
        fd  = &
             physics_get_flow_direction(phys,stat_info_sub)
        fw  = &
             physics_get_flow_width(phys,stat_info_sub)
        
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
        ! Get position, velocity, mass  and  species ID 
        ! of real particles.
        !----------------------------------------------------
        
        CALL particles_get_x(d_particles,x,num_part,stat_info_sub)
        CALL particles_get_v(d_particles,v,num_part,stat_info_sub)
        CALL particles_get_m(d_particles,m,num_part,stat_info_sub)
        CALL particles_get_sid(d_particles,sid,num_part,stat_info_sub)
        
        !----------------------------------------------------
        ! Loop over particles on this process.
        ! Calculate average velocity : V = sum(vi*mi)/N on
        ! this process.
        !----------------------------------------------------
        
        DO j = 1, num_part
           
           !-------------------------------------------------
           ! Count only fluid particles.
           !-------------------------------------------------
           
           IF( sid(j) == 0 ) THEN
              
              !----------------------------------------------
              ! Count only fluid particles in certain width.
              !----------------------------------------------
              
              IF( x(fd,j) < min_phys(fd) + fw .OR. &
                   x(fd,j) > max_phys(fd) - fw ) THEN
                 
                 !-------------------------------------------
                 ! Count mass particles which are involved 
                 ! for calculating flow rate on this process.
                 !-------------------------------------------
                 
                 v_mass  = v_mass + m(j)
                 
                 !-------------------------------------------
                 ! Add involved particles' velocity*mass.
                 !-------------------------------------------
                 
                 this%v_aver(1:num_dim) = &
                      this%v_aver(1:num_dim) + &
                      m(j) * v(1:num_dim,j)
                 
              END IF ! x(fd,j) < fw
              
           END IF ! sid(j) == 0
           
        END DO ! j = 1, num_part_real
        
        
        
        
#ifdef __MPI
        
        !----------------------------------------------------
        ! Calculation in the context of MPI.
        !----------------------------------------------------
        
        comm     = technique_get_comm(tech,stat_info_sub)
        MPI_PREC = technique_get_MPI_PREC(tech,stat_info_sub)
        
        !----------------------------------------------------
        ! Sum up total flow v from all the processes,
        ! then broadcast to all processes.
        !----------------------------------------------------
        
        CALL MPI_ALLREDUCE (this%v_aver(1:num_dim),   &
             v_aver_tot(1:num_dim),num_dim, MPI_PREC, &
             MPI_SUM,comm,stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "statistic_compute_average_v : ", &
                "MPI_ALLREDUCE() for flow rate has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Sum up involved mass of particles for flow rate
        ! from all the processes, 
        ! then broadcast to all processes.
        !----------------------------------------------------
        
        CALL MPI_ALLREDUCE (v_mass,v_mass_tot,1,MPI_PREC, &
             MPI_SUM,comm,stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "statistic_compute_average_v : ", &
                "MPI_ALLREDUCE() for mass  has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! All the processes have the total flow rate.
        !----------------------------------------------------
        
        this%v_aver(1:num_dim)  = v_aver_tot(1:num_dim)
        
        !----------------------------------------------------
        ! All processes have the number of involved particles.
        !----------------------------------------------------
        
        v_mass = v_mass_tot
        
#endif
        
        !----------------------------------------------------
        ! Calculate the average flow rate.
        !----------------------------------------------------
        
        IF ( v_mass /=0.0_MK ) THEN
           
           this%v_aver(1:num_dim) = &
                this%v_aver(1:num_dim) / v_mass
           
        END IF
        
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Release the memory pointed by potiners.
        !----------------------------------------------------

        IF(ASSOCIATED(max_phys)) THEN
           DEALLOCATE(max_phys)
        END IF
        
        IF(ASSOCIATED(min_phys)) THEN
           DEALLOCATE(min_phys)
        END IF
        
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
        
        
        RETURN          
        
      END SUBROUTINE statistic_compute_v_average
      
      
