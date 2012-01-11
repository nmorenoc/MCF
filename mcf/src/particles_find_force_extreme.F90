      SUBROUTINE particles_find_force_extreme(this, &
           comm,MPI_PREC,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_find_force_extreme
        !----------------------------------------------------
        !
        ! Purpose     : Find the minimal and maximal
        !               force per unit mass.
        !               
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     :
        !
        ! Revisions   : V0.1 13.10.2010, original version.
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
        INTEGER, INTENT(IN)                     :: comm
        INTEGER, INTENT(IN)                     :: MPI_PREC
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim, i, num
        REAL(MK), DIMENSION(:), POINTER         :: fa
        REAL(MK)                                :: fa_min
        REAL(MK)                                :: fa_max
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        

        dim = this%num_dim
        
        NULLIFY(fa)
        
        ALLOCATE(fa(this%num_part_real))
        
        num = 0
        
        DO i = 1, this%num_part_real
           
           IF ( this%id(this%sid_idx, i) == &
                mcf_particle_type_fluid ) THEN
              
              num = num+1
              fa(num) = &
                   SQRT(DOT_PRODUCT(this%f(1:dim,i), this%f(1:dim,i)))
              
           END IF
           
        END DO
        
        this%fa_min = &
             MINVAL(fa(1:num))
        this%fa_max = &
             MAXVAL(fa(1:num))
        
        
#ifdef __MPI
        
        !----------------------------------------------------
        ! In context of MPI, pick up the minimum and
        ! maximum from all processes,
        ! then broadcast to all processes.
        !----------------------------------------------------

        CALL MPI_ALLREDUCE (this%fa_min,fa_min, &
             1,MPI_PREC,MPI_MIN,comm,stat_info_sub)
        
        CALL MPI_ALLREDUCE (this%fa_max,fa_max, &
             1,MPI_PREC,MPI_MAX,comm,stat_info_sub)
        
        this%fa_min = fa_min
        this%fa_max = fa_max
        
#endif
        
        IF (ASSOCIATED(fa) ) THEN
           DEALLOCATE(fa)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_find_force_extreme
      
      
      
