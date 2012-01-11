      SUBROUTINE boundary_collect_colloid_interaction(this,&
           comm,MPI_PREC,drag,stat_info)
        !----------------------------------------------------
        ! Subroutine  : boundary_collect_colloid_interaction
        !----------------------------------------------------
        !
        ! Purpose     : Collect total drag/force exerted
        !               on boundaries of colloid
        !               from all processes.
        !                
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     :
        !
        ! Revisions   : V0.1, 19.11.2010, original version.
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
        ! Arguments :
        !----------------------------------------------------
        
        TYPE(Boundary), INTENT(OUT)             :: this
        INTEGER, INTENT(IN)                     :: comm
        INTEGER, INTENT(IN)                     :: MPI_PREC
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        INTEGER                                 :: num
        REAL(MK), DIMENSION(:,:), POINTER       :: t_drag
     
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info        = 0
        stat_info_sub    = 0
        
        dim  = this%num_dim
        num  = 2 * dim

        IF ( SIZE(drag,1) /= dim ) THEN
           PRINT *, "boundary_collect_colloid_interaction : ", &
                "drag dimension does not match ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( SIZE(drag,2) /= num ) THEN
           PRINT *, "boundary_collect_colloid_interaction : ", &
                "drag number does not match ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        NULLIFY(t_drag)
        
#ifdef __MPI
        
        !----------------------------------------------------
        ! In context of MPI, sum up from all processes,
        ! then broadcast to all processes.
        !----------------------------------------------------
        
        ALLOCATE(t_drag(dim,num))
        t_drag(:,:) = 0.0_MK
        
        CALL MPI_ALLREDUCE (drag(:,:), t_drag(:,:),&
             SIZE(drag),MPI_PREC,MPI_SUM,comm,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "boundary_collect_particles_interaction : ", &
                "MPI_ALLREDUCE() for total drag has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%drag(1:dim,1:num) = &
             this%drag(1:dim,1:num) + t_drag(1:dim,1:num)
        
#else
        
        this%drag(1:dim,1:num) = &
             this%drag(1:dim,1:num) + drag(1:dim,1:num)
        
#endif
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(t_drag)) THEN
           DEALLOCATE(t_drag)
        END IF
        
        RETURN          
        
      END SUBROUTINE boundary_collect_colloid_interaction
      
      
