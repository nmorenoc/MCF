      SUBROUTINE io_write_statistic(this,&
           rank,step,time,d_statistic,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_write_statistic
        !----------------------------------------------------
        !
        ! Purpose     : Write statistics information into 
        !               output files.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.2 10.12 2009, change to output
        !               average flow velocity only when
        !               we have fixed_flow_v required from
        !               control component.
        !
        !               V0.1 13.01 2009, original version.
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
        LOGICAL                         :: flow_v_fixed
        LOGICAL                         :: l_p_energy
        INTEGER                         :: num_dim
        REAL(MK)                        :: k_energy
        REAL(MK), DIMENSION(:), POINTER :: momentum
        REAL(MK), DIMENSION(:), POINTER :: v_aver
#ifdef __IO_STATISTIC_STRESS
        REAL(MK), DIMENSION(:), POINTER :: stress
#endif
        REAL(MK)                        :: p_energy
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
        NULLIFY(v_aver)
#ifdef __IO_STATISTIC_STRESS
        NULLIFY(stress)
#endif

        IF ( rank /= 0 ) THEN
           PRINT *, "io_write_statistic : ", &
                "can only be used by root processor !"
           stat_info = -1
           GOTO 9999
        END IF
        
        dynamic_density_ref = &
             control_get_dynamic_density_ref(this%ctrl,stat_info_sub)
        flow_v_fixed = &
             control_get_flow_v_fixed(this%ctrl,stat_info_sub)
        l_p_energy   = &
             control_get_p_energy(this%ctrl,stat_info_sub)
        
        num_dim      = &
             statistic_get_num_dim(d_statistic,stat_info_sub)
        CALL  statistic_get_statistic(d_statistic, &
             k_energy, momentum,stat_info_sub)
        
#ifdef __IO_STATISTIC_STRESS
        CALL statistic_get_stress(d_statistic, & stress, stat_info_sub)
#endif
        !----------------------------------------------------
        ! Kienetic energy and momentum.
        !----------------------------------------------------
        
        data(1)             = k_energy
        data(2:2+num_dim-1) = momentum(1:num_dim)
        
        num_data = 1+num_dim
        
        IF ( flow_v_fixed ) THEN
           
           CALL  statistic_get_v_average(d_statistic, &
                v_aver,stat_info_sub)
           
           data(num_data+1:num_data+num_dim) = &
                v_aver(1:num_dim)
           num_data = num_data + num_dim
           
        END IF
        
        IF( dynamic_density_ref ) THEN
           
           data(num_data+1) = &
                statistic_get_rho_min(d_statistic,stat_info_sub)
           num_data = num_data +1
           
           data(num_data+1) = &
                statistic_get_rho_max(d_statistic,stat_info_sub)
           num_data = num_data +1
           
        END IF
        
#ifdef __IO_STATISTIC_STRESS
        !----------------------------------------------------
        ! For 2D case, off diagonal components only.
        !----------------------------------------------------
        data(num_data+1) = stress(2)
        num_data = num_data + 1
        data(num_data+1) = stress(3)
        num_data = num_data + 1
#endif
        
        IF( l_p_energy ) THEN
           
           p_energy = &
                statistic_get_p_energy(d_statistic,stat_info_sub)
           
           data(num_data+1) = p_energy
           num_data = num_data +1
           
        END IF
        
        WRITE(form, '(A4,I1,A6)') '(I9,', num_data+1, 'E16.8)'
        iform = LEN_TRIM(form)
        
        WRITE(cbuf,form(1:iform)), step, time, data(1:num_data)
        

        WRITE(UNIT=this%statistic_unit,FMT='(A)', &
             IOSTAT=stat_info_sub) TRIM(cbuf)

        IF(stat_info_sub /= 0) THEN
           PRINT *,"io_write_statistic : ",&
                "Writting into statis file failed!"
           stat_info = -1
           GOTO 9999
	END IF
 
 
9999 CONTINUE
        
        IF(ASSOCIATED(momentum)) THEN
           DEALLOCATE(momentum)
        END IF
        
        IF(ASSOCIATED(v_aver)) THEN
           DEALLOCATE(v_aver)
        END IF

#ifdef __IO_STATISTIC_STRESS
        IF(ASSOCIATED(stress)) THEN
           DEALLOCATE(stress)
        END IF
#endif

        RETURN
        
      END SUBROUTINE io_write_statistic
      
