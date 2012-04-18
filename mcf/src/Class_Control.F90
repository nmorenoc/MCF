      MODULE Class_Control
        !----------------------------------------------------
      	! Class	      :	Control
	!----------------------------------------------------
      	!
      	!  Purpose    :
	!> \brief       Variables and corresponding operations
        !>              for control infomation.
      	!>	   	
        !>              The variable memebers are public 
	! Remarks    :
      	!
      	! References :
     	!
      	! Revisions  : 0.2 04.03.2010, including relax_run. 
        !
        !              0.1 03.03.2009, original version.
	!----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        USE mcf_header
        
        IMPLICIT NONE
        SAVE
        
        TYPE Control
           PRIVATE
           
           CHARACTER(LEN=MAX_CHAR)    :: job_name
           CHARACTER(LEN=MAX_CHAR)    :: job_submit_date
           CHARACTER(LEN=10)          :: job_execute_date
           CHARACTER(LEN=10)          :: job_execute_time
           CHARACTER(LEN=10)          :: job_execute_zone
           INTEGER                    :: debug_flag
           LOGICAL                    :: relax_run
           LOGICAL                    :: read_external
           INTEGER                    :: kernel_type
           LOGICAL                    :: symmetry
           INTEGER                    :: rhs_density_type
           LOGICAL                    :: dynamic_density_ref
           INTEGER                    :: stateEquation_type
           LOGICAL                    :: Newtonian
           LOGICAL                    :: Brownian
           INTEGER                    :: random_seed
           INTEGER                    :: rhs_force_type
           LOGICAL                    :: pp_interact_cc
           LOGICAL                    :: pp_interact_cw
           INTEGER                    :: cc_lub_type
           INTEGER                    :: cc_repul_type
           INTEGER                    :: cw_lub_type
           INTEGER                    :: cw_repul_type
           LOGICAL                    :: stress_tensor
           LOGICAL                    :: p_energy
           LOGICAL                    :: flow_v_fixed
           INTEGER                    :: integrate_type
           INTEGER                    :: integrate_colloid_type
           INTEGER                    :: adaptive_dt
           INTEGER                    :: write_output
           INTEGER                    :: write_restart
           
        END TYPE Control
        
        
        INTERFACE control_new
           MODULE PROCEDURE control_init_default
        END INTERFACE
        
      CONTAINS
        
#include "control_new.F90"
#include "control_finalize.F90"
#include "control_check_parameters.F90"
#include "control_get.F90"
#include "control_set.F90"
        
          
      END MODULE Class_Control
      
