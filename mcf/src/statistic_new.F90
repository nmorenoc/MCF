      SUBROUTINE statistic_init(this,d_ctrl,d_num_dim,stat_info)

        TYPE(Statistic), INTENT(OUT)          :: this
        TYPE(Control), INTENT(IN)             :: d_ctrl
        INTEGER, INTENT(IN)                   :: d_num_dim
        INTEGER, INTENT(OUT)                  :: stat_info
        
        
        INTEGER                               :: stat_info_sub
        
        
        stat_info     = 0
        stat_info_sub = 0
        
        this%dynamic_density_ref = &
             control_get_dynamic_density_ref(d_ctrl,stat_info_sub)
        this%flow_v_fixed = &
             control_get_flow_v_fixed(d_ctrl,stat_info_sub)
        this%l_p_energy   = &
             control_get_p_energy(d_ctrl,stat_info_sub)
        
        IF ( d_num_dim < 2 .OR. &
             d_num_dim > 3 ) THEN
           PRINT *, "statistic_init : ", &
                "Dimension not supported !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%num_dim    = d_num_dim
        
        this%k_energy              = 0.0_MK
        this%momentum(1:d_num_dim) = 0.0_MK
        this%disorder(1:d_num_dim) = 0.0_MK
        this%v_aver(1:d_num_dim)   = 0.0_MK
        this%stress(1:d_num_dim**2)= 0.0_MK
        this%p_energy              = 0.0_MK
        this%rho_min               = 0.0_MK
        this%rho_max               = 0.0_MK

        this%colloid_implicit_pair_sweep_adaptive = .FALSE.
        this%colloid_implicit_pair_num_sweep      = 1
        this%colloid_implicit_pair_sweep_error    = 1.0e2_MK
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE statistic_init
      
      
      SUBROUTINE statistic_display_parameters(this,stat_info)
        
        TYPE(Statistic), INTENT(IN)     :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        PRINT *,'--------------------------------------'
        PRINT *,'      Statistic Parameters '
        PRINT *,'--------------------------------------'
        
        PRINT *, "num_dim          : ", this%num_dim
        PRINT *, "Kinetic energy   : ", this%k_energy
        PRINT *, "Totoal momentum  : ", this%momentum(1:this%num_dim)
        PRINT *, "Disorder level   : ", this%disorder(1:this%num_dim)
        PRINT *, "Average velocity : ", this%v_aver(1:this%num_dim)
        PRINT *, "Totoal stress    : ", this%stress(1:this%num_dim**2)
        
        IF( this%l_p_energy ) THEN
           PRINT *, "Potential energy : ", "Needed "
           PRINT *, "p_energy         : ", this%p_energy
        ELSE
           PRINT *, "Potential energy : ", "Not Needed "
        END IF
        
        RETURN          
        
      END SUBROUTINE statistic_display_parameters
      
