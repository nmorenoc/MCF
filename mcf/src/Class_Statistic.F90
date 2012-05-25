      MODULE Class_Statistic
        !----------------------------------------------------
      	!  Class      :	Statistic
	!----------------------------------------------------
      	!
      	!  Purpose    :
	!> \brief       Variables and corresponding operations
        !>              for statistic quantities.
      	!>	   	
        !>              The variable memebers are public 
        !
        !  Remarks    : Since there is a 'SAVE' after 'IMPLICIT NONE',
        !               all procedures using this Module sharing the
        !               same copy of this Module.
        !
      	!
      	!  References :
     	!
      	!  Revisions  : 0.1 03.03.2009
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
             
        USE Class_Control
        USE Class_Technique
        USE Class_Particles
        USE Class_Debug
        
        IMPLICIT NONE
        SAVE
        
        TYPE Statistic

           PRIVATE
           TYPE(Control)                      :: ctrl
           LOGICAL                            :: dynamic_density_ref
           LOGICAL                            :: flow_v_fixed
           LOGICAL                            :: l_p_energy
           INTEGER                            :: num_dim
           REAL(MK)                           :: k_energy
           REAL(MK), DIMENSION(3)             :: momentum
           REAL(MK), DIMENSION(3)             :: disorder
           REAL(MK), DIMENSION(3)             :: v_aver
           REAL(MK)                           :: rho_min
           REAL(MK)                           :: rho_max
           REAL(MK), DIMENSION(9)             :: stress
           REAL(MK), DIMENSION(9)             :: stress_p
           REAL(MK), DIMENSION(9)             :: stress_v
           REAL(MK), DIMENSION(9)             :: stress_r
           REAL(MK)                           :: p_energy
           LOGICAL                            :: colloid_implicit_pair_sweep_adaptive
           INTEGER                            :: colloid_implicit_pair_num_sweep
           REAL(MK)                           :: colloid_implicit_pair_sweep_error
           
        END TYPE Statistic
        
        INTERFACE statistic_new
           MODULE PROCEDURE statistic_init
        END INTERFACE
        
        CONTAINS
          
#include "statistic_new.F90"
#include "statistic_finalize.F90"
#include "statistic_compute_disorder.F90"
#include "statistic_compute_statistic.F90"
#include "statistic_compute_v_average.F90"
#include "statistic_get.F90"
#include "statistic_set.F90"

        END MODULE Class_Statistic
