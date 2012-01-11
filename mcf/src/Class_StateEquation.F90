      MODULE Class_StateEquation
        !----------------------------------------------------
      	!  Class      :	StateEquation
	!----------------------------------------------------
      	!
      	!  Purpose    :
	!> \brief       Variables and corresponding operations
        !>              for state equation quantities.
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
           
        USE mcf_header
        
        IMPLICIT NONE
        SAVE
        
        TYPE StateEquation
           INTEGER              :: stateEquation_type
           REAL(MK)             :: c 
           REAL(MK)             :: p0
           REAL(MK)             :: rho_ref
           REAL(MK)             :: gamma

        END TYPE StateEquation
        
        INTERFACE stateEquation_new
           MODULE PROCEDURE stateEquation_init
        END INTERFACE
        
        INTERFACE stateEquation_compute_pressure
           MODULE PROCEDURE stateEquation_compute_pressure_scalar
           MODULE PROCEDURE stateEquation_compute_pressure_vector
        END INTERFACE
        
      CONTAINS
          
#include "stateEquation_new.F90"
#include "stateEquation_finalize.F90"
#include "stateEquation_compute_pressure.F90"
#include "stateEquation_get.F90"
#include "stateEquation_set.F90"

      END MODULE Class_StateEquation
