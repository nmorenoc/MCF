      SUBROUTINE boundary_collect_particles_interaction(this,&
           comm,MPI_PREC,drag,drag_p,drag_v,drag_r,stat_info)
        !----------------------------------------------------
        ! Subroutine  : boundary_collect_particles_interaction
        !----------------------------------------------------
        !
        ! Purpose     : Collect total drag/force exerted
        !               on boundaries from all processes.
        !                
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     :
        !
        ! Revisions   : V0.1, 23.10.2009, original version.
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
        REAL(MK),DIMENSION(:,:),INTENT(IN)      :: drag
        REAL(MK),DIMENSION(:,:),INTENT(IN),OPTIONAL:: drag_p
        REAL(MK),DIMENSION(:,:),INTENT(IN),OPTIONAL:: drag_v
        REAL(MK),DIMENSION(:,:),INTENT(IN),OPTIONAL:: drag_r
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        INTEGER                                 :: num
        REAL(MK), DIMENSION(:,:), POINTER       :: t_drag
#ifdef __FORCE_SEPARATE
        REAL(MK), DIMENSION(:,:), POINTER       :: t_drag_p
        REAL(MK), DIMENSION(:,:), POINTER       :: t_drag_v
        REAL(MK), DIMENSION(:,:), POINTER       :: t_drag_r
#endif
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info        = 0
        stat_info_sub    = 0
        
        dim  = this%num_dim
        num  = 2 * dim

        IF ( SIZE(drag,1) /= dim ) THEN
           PRINT *, "boundary_collect_particles_interaction : ", &
                "drag dimension does not match ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( SIZE(drag,2) /= num ) THEN
           PRINT *, "boundary_collect_particles_interaction : ", &
                "drag number does not match ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        this%drag(:,:) = 0.0_MK        
        NULLIFY(t_drag)

#ifdef __FORCE_SEPARATE
        this%drag_p(:,:) = 0.0_MK        
        NULLIFY(t_drag_p)
        this%drag_v(:,:) = 0.0_MK        
        NULLIFY(t_drag_v)
        this%drag_r(:,:) = 0.0_MK        
        NULLIFY(t_drag_r)
#endif
        
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
        
        this%drag(1:dim,1:num) = t_drag(1:dim,1:num)
        
#ifdef __FORCE_SEPARATE

        ALLOCATE(t_drag_p(dim,num))
        t_drag_p(:,:) = 0.0_MK        
        CALL MPI_ALLREDUCE (drag_p(:,:), t_drag_p(:,:),&
             SIZE(drag_p),MPI_PREC,MPI_SUM,comm,stat_info_sub)        
        this%drag_p(1:dim,1:num) = t_drag_p(1:dim,1:num)

        ALLOCATE(t_drag_v(dim,num))
        t_drag_v(:,:) = 0.0_MK        
        CALL MPI_ALLREDUCE (drag_v(:,:), t_drag_v(:,:),&
             SIZE(drag_v),MPI_PREC,MPI_SUM,comm,stat_info_sub)        
        this%drag_v(1:dim,1:num) = t_drag_v(1:dim,1:num)
        
        ALLOCATE(t_drag_r(dim,num))
        t_drag_r(:,:) = 0.0_MK        
        CALL MPI_ALLREDUCE (drag_r(:,:), t_drag_r(:,:),&
             SIZE(drag_r),MPI_PREC,MPI_SUM,comm,stat_info_sub)        
        this%drag_r(1:dim,1:num) = t_drag_r(1:dim,1:num)
        
#endif

#else
        
        this%drag(1:dim,1:num) = drag(1:dim,1:num)
        
#ifdef __FORCE_SEPARATE
        
        this%drag_p(1:dim,1:num) = drag_p(1:dim,1:num)
        this%drag_v(1:dim,1:num) = drag_v(1:dim,1:num)
        this%drag_r(1:dim,1:num) = drag_r(1:dim,1:num)
        
#endif

#endif
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(t_drag)) THEN
           DEALLOCATE(t_drag)
        END IF
        
#ifdef __FORCE_SEPARATE
        IF(ASSOCIATED(t_drag_p)) THEN
           DEALLOCATE(t_drag_p)
        END IF
        IF(ASSOCIATED(t_drag_v)) THEN
           DEALLOCATE(t_drag_v)
        END IF
        IF(ASSOCIATED(t_drag_r)) THEN
           DEALLOCATE(t_drag_r)
        END IF
#endif
        
        RETURN          
        
      END SUBROUTINE boundary_collect_particles_interaction
      
      
