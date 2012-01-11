      SUBROUTINE io_write_statistic_relax(this,&
           rank,step,time,d_statistic,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_write_statistic_relax
        !----------------------------------------------------
        !
        ! Purpose     : Write statistics and relax information
        !               into output files.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 09.03 2010, original version.
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
	! Arguments :        
        !----------------------------------------------------
        
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)             :: rank
        INTEGER,  INTENT(IN)	        :: step
        REAL(MK), INTENT(IN)	        :: time
        TYPE(Statistic), INTENT(IN)     :: d_statistic
        INTEGER,  INTENT(OUT)	        :: stat_info
        
        !----------------------------------------------------
	! Local variables start here :        
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: dynamic_density_ref
        INTEGER                         :: num_dim
        REAL(MK)                        :: k_energy
        REAL(MK), DIMENSION(:), POINTER :: momentum
        REAL(MK), DIMENSION(:), POINTER :: disorder
        INTEGER                         :: num_data
        REAL(MK), DIMENSION(8)          :: data
        CHARACTER(len=MAX_CHAR)	        :: form
        INTEGER                         :: iform
        CHARACTER(len=MAX_CHAR)	        :: cbuf
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(momentum)
        NULLIFY(disorder)
        
        IF ( rank /= 0 ) THEN
           PRINT *, "io_write_statistic_relax : ", &
                "can only be used by root processor !"
           stat_info = -1
           GOTO 9999
        END IF
        
        dynamic_density_ref = &
             control_get_dynamic_density_ref(this%ctrl,stat_info_sub)
        num_dim      = &
             statistic_get_num_dim(d_statistic,stat_info_sub)
        CALL  statistic_get_statistic(d_statistic, &
             k_energy, momentum,stat_info_sub)
        
        !----------------------------------------------------
        ! Kienetic energy and momentum.
        !----------------------------------------------------
        
        data(1)             = k_energy
        data(2:2+num_dim-1) = momentum(1:num_dim)
        
        num_data = 1+num_dim
        
        CALL  statistic_get_disorder(d_statistic, &
             disorder,stat_info_sub)
        
        data(num_data+1:num_data+num_dim) = &
             disorder(1:num_dim)
        
        num_data = num_data + num_dim
        
        IF( dynamic_density_ref ) THEN
           
           data(num_data+1) = &
                statistic_get_rho_min(d_statistic,stat_info_sub)
           num_data = num_data +1
           
           data(num_data+1) = &
                statistic_get_rho_max(d_statistic,stat_info_sub)
           num_data = num_data +1
           
        END IF
        
        WRITE(form, '(A4,I1,A6)') '(I9,', num_data+1, 'E16.8)'
        iform = LEN_TRIM(form)
        
        WRITE(cbuf,form(1:iform)), step, time, data(1:num_data)

        
        WRITE(UNIT=this%statistic_relax_unit,FMT='(A)', &
             IOSTAT=stat_info_sub) TRIM(cbuf)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *,"io_write_statistic_relax : ",&
                "Writting into statis file failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
9999    CONTINUE
        
        IF(ASSOCIATED(momentum)) THEN
           DEALLOCATE(momentum)
        END IF
        
        IF(ASSOCIATED(disorder)) THEN
           DEALLOCATE(disorder)
        END IF

        RETURN
        
      END SUBROUTINE io_write_statistic_relax
