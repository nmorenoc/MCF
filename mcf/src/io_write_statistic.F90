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
        LOGICAL                         :: Brownian
        LOGICAL                         :: dynamic_density_ref
        LOGICAL                         :: flow_v_fixed
        LOGICAL                         :: l_stress_tensor
        LOGICAL                         :: l_p_energy
        LOGICAL                         :: integrate_colloid_type
        
        INTEGER                         :: num_dim
        REAL(MK)                        :: k_energy
        REAL(MK), DIMENSION(:), POINTER :: momentum
        REAL(MK), DIMENSION(:), POINTER :: v_aver
        REAL(MK), DIMENSION(:), POINTER :: stress
#ifdef __IO_STATISTIC_STRESS_SEPARATE
        REAL(MK), DIMENSION(:), POINTER :: stress_p
        REAL(MK), DIMENSION(:), POINTER :: stress_v
        REAL(MK), DIMENSION(:), POINTER :: stress_r
#endif
        REAL(MK)                        :: p_energy

        LOGICAL                         :: coll_implicit_pair_sweep_adaptive
        INTEGER                         :: coll_implicit_pair_num_sweep
        REAL(MK)                        :: coll_implicit_pair_sweep_error
     
        INTEGER                         :: num_data
        REAL(MK), DIMENSION(20)         :: data
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
        NULLIFY(stress)

#ifdef __IO_STATISTIC_STRESS_SEPARATE
        NULLIFY(stress_p)
        NULLIFY(stress_v)
        NULLIFY(stress_r)
#endif
        
        IF ( rank /= 0 ) THEN
           PRINT *, "io_write_statistic : ", &
                "can only be used by root processor !"
           stat_info = -1
           GOTO 9999
        END IF

        Brownian = &
             control_get_Brownian(this%ctrl,stat_info_sub)     
        dynamic_density_ref = &
             control_get_dynamic_density_ref(this%ctrl,stat_info_sub)
        flow_v_fixed = &
             control_get_flow_v_fixed(this%ctrl,stat_info_sub)
        l_stress_tensor = &
             control_get_stress_tensor(this%ctrl,stat_info_sub)
        l_p_energy   = &
             control_get_p_energy(this%ctrl,stat_info_sub)
        integrate_colloid_type = &
             control_get_integrate_colloid_type(this%ctrl,stat_info_sub)
        
        num_dim      = &
             statistic_get_num_dim(d_statistic,stat_info_sub)
        CALL  statistic_get_statistic(d_statistic, &
             k_energy, momentum,stat_info_sub)
        
        IF ( l_stress_tensor ) THEN
           
           CALL statistic_get_stress(d_statistic, stress, stat_info_sub)
           
#ifdef __IO_STATISTIC_STRESS_SEPARATE
           CALL statistic_get_stress_p(d_statistic, stress_p, stat_info_sub)
           CALL statistic_get_stress_v(d_statistic, stress_v, stat_info_sub)
           IF ( Brownian ) THEN
              CALL statistic_get_stress_r(d_statistic, stress_r, stat_info_sub)
           END IF
#endif
           
        END IF
        
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

        else

#ifdef __IO_STATISTIC_V_AVERAGE
           
           CALL  statistic_get_v_average(d_statistic, &
                v_aver,stat_info_sub)
           
           data(num_data+1:num_data+num_dim) = &
                v_aver(1:num_dim)
           num_data = num_data + num_dim           
#endif
           
        END IF
        
        IF( dynamic_density_ref ) THEN
           
           data(num_data+1) = &
                statistic_get_rho_min(d_statistic,stat_info_sub)
           num_data = num_data + 1
           
           data(num_data+1) = &
                statistic_get_rho_max(d_statistic,stat_info_sub)
           num_data = num_data + 1
           
        END IF
        
        !----------------------------------------------------
        ! For 2D case, off diagonal components only.
        !----------------------------------------------------
        
        IF ( l_stress_tensor ) THEN
           
           data(num_data+1) = stress(2)
           num_data = num_data + 1
           data(num_data+1) = stress(3)
           num_data = num_data + 1
           
#ifdef __IO_STATISTIC_STRESS_SEPARATE           
           data(num_data+1) = stress_p(2)
           num_data = num_data + 1
           data(num_data+1) = stress_p(3)
           num_data = num_data + 1
           data(num_data+1) = stress_v(2)
           num_data = num_data + 1
           data(num_data+1) = stress_v(3)
           num_data = num_data + 1
           
           IF ( Brownian ) THEN
              data(num_data+1) = stress_r(2)
              num_data = num_data + 1
              data(num_data+1) = stress_r(3)
              num_data = num_data + 1
           END IF           
#endif
           
        END IF
        
        IF( l_p_energy ) THEN
           
           p_energy = &
                statistic_get_p_energy(d_statistic,stat_info_sub)
           
           data(num_data+1) = p_energy
           num_data = num_data + 1
           
        END IF

        
        IF ( integrate_colloid_type == -2 ) THEN
           
           coll_implicit_pair_sweep_adaptive = &
                statistic_get_colloid_implicit_pair_sweep_adaptive(d_statistic,stat_info_sub)
           
           IF ( coll_implicit_pair_sweep_adaptive ) THEN

              coll_implicit_pair_num_sweep = &
                   statistic_get_colloid_implicit_pair_num_sweep(d_statistic,stat_info_sub)
              
              coll_implicit_pair_sweep_error = &
                   statistic_get_colloid_implicit_pair_sweep_error(d_statistic,stat_info_sub)
              
              data(num_data+1) = coll_implicit_pair_num_sweep
              data(num_data+2) = coll_implicit_pair_sweep_error
         
              num_data = num_data + 2
              
           END IF

        END IF
        
        !----------------------------------------------------
        ! For the format, it depends on number of columns.
        !----------------------------------------------------
        IF ( num_data + 1 < 10 ) THEN
           WRITE(form, '(A4,I1,A6)') '(I9,', num_data+1, 'E16.8)'
        ELSE
           WRITE(form, '(A4,I2,A6)') '(I9,', num_data+1, 'E16.8)'
        END IF
        
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
        
        IF(ASSOCIATED(stress)) THEN
           DEALLOCATE(stress)
        END IF

#ifdef __IO_STATISTIC_STRESS_SEPARATE
        IF(ASSOCIATED(stress_p)) THEN
           DEALLOCATE(stress_p)
        END IF
        IF(ASSOCIATED(stress_v)) THEN
           DEALLOCATE(stress_v)
        END IF
        IF(ASSOCIATED(stress_r)) THEN
           DEALLOCATE(stress_r)
        END IF        
#endif


        
        RETURN
        
      END SUBROUTINE io_write_statistic
      
