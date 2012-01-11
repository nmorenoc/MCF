      SUBROUTINE statistic_compute_disorder(this,d_particles,&
           dx0,stat_info)
        !----------------------------------------------------
        ! Subroutine  : statistic_compute_disorder
        !----------------------------------------------------
        !
        ! Purpose     : Computes disorder level of particles.
        !               
        !
        ! Routines    :
        !
        ! Remarks     :
        !
        ! References  :
        !
        ! Revisions   : V0.1 09.03 2010, original version.
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
        ! dx0           :  initial particles' distances.
        ! stat_info     :  status information.
        !----------------------------------------------------
        
        TYPE(Statistic), INTENT(INOUT)          :: this
        TYPE(Particles), INTENT(IN)             :: d_particles
        REAL(MK), DIMENSION(:), INTENT(IN)      :: dx0
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num_dim

        !----------------------------------------------------
        ! Particles variables,  
        ! position, velocity, mass and species ID.
        !----------------------------------------------------

        INTEGER                                 :: num_part
        REAL(MK), DIMENSION(:,:), POINTER       :: x
        INTEGER,  DIMENSION(:), POINTER         :: sid

        !----------------------------------------------------
        ! A technique object pointer.
        !----------------------------------------------------
        
        TYPE(Technique), POINTER                :: tech
        INTEGER                                 :: comm
        INTEGER                                 :: MPI_PREC
        
        !----------------------------------------------------
        ! Total values of disorder.
        !----------------------------------------------------
        
        REAL(MK), DIMENSION(3)                  :: disorder_tot
        INTEGER                                 :: n, n_tot
        
        !----------------------------------------------------
        ! k    : coefficient
        ! i,j  : Integer counters and indices.
        !----------------------------------------------------
        
        REAL(MK),DIMENSION(3)                   :: k
        INTEGER                                 :: j
        

        !----------------------------------------------------
        ! Check if dimensions match.
        !----------------------------------------------------
        
        num_dim = this%num_dim
        
        IF(SIZE(dx0,1) /= this%num_dim) THEN
           PRINT *, "statistic_compute_disorder : ", &
                "dx0 has different dimension !"
           stat_info = -1
           GOTO 9999           
        END IF
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(x)
        NULLIFY(sid)
        
        NULLIFY(tech)
        
           
        k(1:num_dim) = &
             2.0_MK * mcf_pi / dx0(1:num_dim)
        
        
        this%disorder(1:num_dim) = 0.0_MK
        disorder_tot(1:num_dim) = 0.0_MK
        
        n     = 0
        n_tot = 0
        
        
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
        ! Get position and sid of particles.
        !----------------------------------------------------
        
        CALL particles_get_x(d_particles,x,num_part,stat_info_sub)
        CALL particles_get_sid(d_particles,sid,num_part,stat_info_sub)

        !----------------------------------------------------
        ! Disorder r -= COS(2*pi*x/dx0)
        !          r /= n;
        !----------------------------------------------------
        
        DO j = 1, num_part
           
           !-------------------------------------------------
           ! Count only fluid particles.
           !-------------------------------------------------
           
           IF( sid(j) == 0 ) THEN
              
              
              this%disorder(1:num_dim) = &
                   this%disorder(1:num_dim) - &
                   COS(k(1:num_dim)*x(1:num_dim,j))
                 
              
              n = n + 1
              
           END IF ! sid(j) == 0
           
        END DO ! j = 1, num_part
        
        !----------------------------------------------------
        ! Calculation in the context of MPI.
        !----------------------------------------------------
        
#ifdef __MPI
        
        comm     = technique_get_comm(tech,stat_info_sub)
        MPI_PREC = technique_get_MPI_PREC(tech,stat_info_sub)
        
        !----------------------------------------------------
        ! Sum up the disorder from all the processes,
        ! broadcast the result.
        !----------------------------------------------------
        
        CALL MPI_ALLREDUCE (this%disorder(1:num_dim),  &
             disorder_tot(1:num_dim),num_dim,MPI_PREC, &
             MPI_SUM,comm,stat_info_sub) 
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "statistic_compute_disorder : ", &
                "MPI_ALLREDUCE() for disorder has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Sum up the number of particles involved from all
        ! processes, broadcast the result.
        !----------------------------------------------------
        
        CALL MPI_ALLREDUCE (n, n_tot,1,MPI_INTEGER, &
             MPI_SUM,comm,stat_info_sub) 
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "statistic_compute_disorder : ", &
                "MPI_ALLREDUCE() for n has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! All processes get total disorder and n.
        !----------------------------------------------------
        
        this%disorder(1:num_dim) = disorder_tot(1:num_dim)
        n   = n_tot
        
#endif
        
        !----------------------------------------------------
        ! Calculate the disorder
        !----------------------------------------------------
        
        IF ( n /=0 ) THEN
           
           this%disorder(1:num_dim) = &
                this%disorder(1:num_dim) / n
           
        END IF

9999    CONTINUE
        
        !----------------------------------------------------
        ! Release the memory pointed by potiners.
        !----------------------------------------------------
        
        IF(ASSOCIATED(x)) THEN
           DEALLOCATE(x)
        END IF
        
        IF(ASSOCIATED(sid)) THEN
           DEALLOCATE(sid)
        END IF
        
        
        RETURN          
        
      END SUBROUTINE statistic_compute_disorder
      
      
