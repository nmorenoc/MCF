!--------------------------------------------------
! Subroutine  :  control_set_*
!--------------------------------------------------
!
! Purpose     : Set routines of Class Control.
!
! Reference   :
!
! Remark      :
!
! Revisions   : V0.1 01.03.2009, original version.
!
!--------------------------------------------------
! Author      : Xin Bian
! Contact     : xin.bian@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!--------------------------------------------------

      SUBROUTINE control_set_job_name(this,d_job_name,stat_info)

        TYPE(Control), INTENT(INOUT)             :: this
        CHARACTER(*)                             :: d_job_name
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info     = 0
        this%job_name = d_job_name
        
        RETURN 
        
      END SUBROUTINE control_set_job_name

      
      SUBROUTINE control_set_job_submit_date(this,d_job_date,stat_info)

        TYPE(Control), INTENT(INOUT)             :: this
        CHARACTER(*)                             :: d_job_date
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info     = 0
        this%job_submit_date = d_job_date
        
        RETURN 
        
      END SUBROUTINE control_set_job_submit_date
      
      
      SUBROUTINE control_set_debug_flag(this,d_debug_flag,stat_info)

        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_debug_flag
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        this%debug_flag =   d_debug_flag 
        
        RETURN 
        
      END SUBROUTINE control_set_debug_flag

      
      SUBROUTINE control_set_relax_run(this,d_relax_run,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        LOGICAL, INTENT(IN)                      :: d_relax_run
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%relax_run =d_relax_run
        
        RETURN
        
      END SUBROUTINE control_set_relax_run


      SUBROUTINE control_set_read_external(this,d_read_external,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        LOGICAL, INTENT(IN)                      :: d_read_external
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%read_external =d_read_external
        
        RETURN
        
      END SUBROUTINE control_set_read_external


      SUBROUTINE control_set_kernel_type(this,d_kernel_type,stat_info)

        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_kernel_type
        INTEGER, INTENT(out)                     :: stat_info

        stat_info = 0
        
        this%kernel_type = d_kernel_type
        
        RETURN

      END SUBROUTINE control_set_kernel_type


      SUBROUTINE control_set_symmetry(this,d_symmetry,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        LOGICAL, INTENT(IN)                      :: d_symmetry
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%symmetry =d_symmetry
        
        RETURN
        
      END SUBROUTINE control_set_symmetry    
   

      SUBROUTINE control_set_rhs_density_type(this,d_rhs_density_type,stat_info)

        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_rhs_density_type
        INTEGER, INTENT(out)                     :: stat_info

        stat_info = 0
        this%rhs_density_type = d_rhs_density_type
        
        RETURN
        
      END SUBROUTINE control_set_rhs_density_type
      

      SUBROUTINE control_set_dynamic_density_ref(this,d_ref,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        LOGICAL, INTENT(IN)                      :: d_ref
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%dynamic_density_ref =d_ref
        
        RETURN
        
      END SUBROUTINE control_set_dynamic_density_ref


      SUBROUTINE control_set_stateEquation_type(this,d_stateE,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_stateE
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        this%stateEquation_type = d_stateE
        
        RETURN
        
      END SUBROUTINE control_set_stateEquation_type
      
      
      SUBROUTINE control_set_Newtonian(this,d_Newtonian,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        LOGICAL, INTENT(IN)                      :: d_Newtonian
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%Newtonian =d_Newtonian
        
        RETURN
        
      END SUBROUTINE control_set_Newtonian

  
      SUBROUTINE control_set_Brownian(this,d_Brownian,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        LOGICAL, INTENT(IN)                      :: d_Brownian
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%Brownian =d_Brownian
        
        RETURN
        
      END SUBROUTINE control_set_Brownian
      

      SUBROUTINE control_set_random_seed(this,d_random_seed,stat_info)

        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_random_seed
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        this%random_seed = d_random_seed
        
        RETURN
        
      END SUBROUTINE control_set_random_seed
     
      
      SUBROUTINE control_set_rhs_force_type(this,d_rhs_force_type,stat_info)

        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_rhs_force_type
        INTEGER, INTENT(out)                     :: stat_info

        stat_info = 0
        this%rhs_force_type = d_rhs_force_type
        
        RETURN
        
      END SUBROUTINE control_set_rhs_force_type
      
      
      SUBROUTINE control_set_pp_interact_cc(this,d_pp_interact_cc,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        LOGICAL, INTENT(IN)                      :: d_pp_interact_cc
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%pp_interact_cc =d_pp_interact_cc
        
        RETURN
        
      END SUBROUTINE control_set_pp_interact_cc
      
  
      SUBROUTINE control_set_pp_interact_cw(this,d_pp_interact_cw,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        LOGICAL, INTENT(IN)                      :: d_pp_interact_cw
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%pp_interact_cw =d_pp_interact_cw
        
        RETURN
        
      END SUBROUTINE control_set_pp_interact_cw
      
      
      SUBROUTINE control_set_cc_lub_type(this,d_lub_type,stat_info)

        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_lub_type
        INTEGER, INTENT(out)                     :: stat_info

        stat_info = 0
        this%cc_lub_type = d_lub_type
        
        RETURN
        
      END SUBROUTINE control_set_cc_lub_type
      
      
      SUBROUTINE control_set_cc_repul_type(this,d_repul_type,stat_info)

        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_repul_type
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        this%cc_repul_type = d_repul_type
        
        RETURN
        
      END SUBROUTINE control_set_cc_repul_type
    
      
      SUBROUTINE control_set_cw_lub_type(this,d_lub_type,stat_info)

        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_lub_type
        INTEGER, INTENT(out)                     :: stat_info

        stat_info = 0
        this%cw_lub_type = d_lub_type
        
        RETURN
        
      END SUBROUTINE control_set_cw_lub_type
      
      
      SUBROUTINE control_set_cw_repul_type(this,d_repul_type,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_repul_type
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        this%cw_repul_type = d_repul_type
        
        RETURN
        
      END SUBROUTINE control_set_cw_repul_type
      

      SUBROUTINE control_set_p_energy(this,d_p_energy,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        LOGICAL, INTENT(IN)                      :: d_p_energy
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%p_energy =d_p_energy
        
        RETURN
        
      END SUBROUTINE control_set_p_energy


      SUBROUTINE control_set_flow_v_fixed(this,d_flow_v_fixed,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        LOGICAL, INTENT(IN)                      :: d_flow_v_fixed
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%flow_v_fixed =d_flow_v_fixed
        
        RETURN
        
      END SUBROUTINE control_set_flow_v_fixed
      
     
      SUBROUTINE control_set_integrate_type(this,d_integrate_type,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_integrate_type
        INTEGER, INTENT(out)                     :: stat_info
       
        stat_info = 0
        
        this%integrate_type = d_integrate_type
       
        RETURN
        
      END SUBROUTINE control_set_integrate_type


      SUBROUTINE control_set_adaptive_dt(this,d_adaptive_dt,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_adaptive_dt
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%adaptive_dt =d_adaptive_dt
        
        RETURN
        
      END SUBROUTINE control_set_adaptive_dt
      

      SUBROUTINE control_set_write_output(this,d_write,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_write
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%write_output =d_write
        
        RETURN
        
      END SUBROUTINE control_set_write_output
      

      SUBROUTINE control_set_write_restart(this,d_write,stat_info)
        
        TYPE(Control), INTENT(INOUT)             :: this
        INTEGER, INTENT(IN)                      :: d_write
        INTEGER, INTENT(out)                     :: stat_info
        
        stat_info = 0
        
        this%write_restart =d_write
        
        RETURN
        
      END SUBROUTINE control_set_write_restart
      
