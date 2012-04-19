      SUBROUTINE colloid_collect_particles_interaction(this,&
           comm,MPI_PREC,drag,torque,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_collect_particles_interaction
        !----------------------------------------------------
        !
        ! Purpose     : Collect force and torque on all 
        !               processes from colloidal boundary
        !               particles. 
        !
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     :
        !
        ! Revisions   :  V0.1 19.11.2010,  original version.
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
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        INTEGER, INTENT(IN)                     :: comm
        INTEGER, INTENT(IN)                     :: MPI_PREC
        REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: drag
        REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: torque
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !
        ! t_drag   : total drag on all processes
        ! t_torque : total torque on all processes        
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim,num
        REAL(MK), ALLOCATABLE, DIMENSION(:,:)   :: t_drag
        REAL(MK), ALLOCATABLE, DIMENSION(:,:)   :: t_torque

        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        dim = this%num_dim
        num = this%num_colloid
        
        IF ( SIZE(drag,1) /= dim ) THEN
           PRINT *, "colloid_collect_particles_interaction: ", &
                "input drag dimension does not match!"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( SIZE(drag,2) /= num ) THEN
           PRINT *, "colloid_collect_particles_interaction: ", &
                "input drag number does not match!"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( SIZE(torque,1) /= 3 ) THEN
           PRINT *, "colloid_collect_particles_interaction: ", &
                "input torque dimension does not match!"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( SIZE(torque,2) /= num ) THEN
           PRINT *, "colloid_collect_particles_interaction: ", &
                "input drag number does not match!"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! collect boundary particles' contribution of
        ! total drag and torque on colloids.
        !----------------------------------------------------

#ifdef __MPI
        
        !----------------------------------------------------
        ! In context of MPI, sum up from all processes,
        ! then broadcast to all processes.
        !
        ! Loop over all colloid particles.
        !----------------------------------------------------
        
        ALLOCATE(t_drag(1:dim,1:num))
        ALLOCATE(t_torque(1:3,1:num))
        
        t_drag(:,:) = 0.0_MK
        t_torque(:,:) = 0.0_MK
        
        CALL MPI_ALLREDUCE (drag(:,:),t_drag(:,:), &
             size(drag),MPI_PREC,MPI_SUM,comm,stat_info_sub)
           
        IF (stat_info_sub /= 0 ) THEN
           PRINT *, "colloid_collect_particles_interaction : ", &
                "MPI_ALLREDUCE() for total drag has problem !"
           stat_info = -1
           GOTO 9999
        END IF
           
        drag(1:dim,1:num) = t_drag(1:dim,1:num)
        
        CALL MPI_ALLREDUCE (torque(:,:),t_torque(:,:), &
             SIZE(t_torque),MPI_PREC,MPI_SUM,comm,stat_info_sub)
        
        IF (stat_info_sub /= 0 ) THEN
           PRINT *, "colloid_collect_particles_interaction : ", &
                "MPI_ALLREDUCE() for total torque has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        torque(1:3,1:num) = t_torque(1:3,1:num)

#else
        
        ! For single processor, nothing needs to be done!
        
#endif
        
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE colloid_collect_particles_interaction
      
      
      
