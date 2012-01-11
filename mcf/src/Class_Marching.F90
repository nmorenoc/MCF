      MODULE Class_Marching
        !----------------------------------------------------
      	!  Class      :	Marching
	!----------------------------------------------------
      	!
      	!  Purpose    :
	!> \brief       Variables and corresponding operations
        !>              for Marching quantities.
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
        USE Class_Debug
        
        USE Class_Control
        USE Class_Physics
        USE Class_Particles
        USE Class_IO
        USE Class_Statistic
        
        IMPLICIT NONE
        SAVE
        
        TYPE Marching
           PRIVATE
           TYPE(Control), POINTER             :: ctrl
           TYPE(Physics), POINTER             :: phys 
           TYPE(IO), POINTER                  :: io
           TYPE(Particles), POINTER           :: particles
           TYPE(Technique), POINTER           :: tech
           TYPE(Statistic)                    :: statis
           INTEGER                            :: integrate_type
        END TYPE Marching
        
        INTERFACE marching_new
           MODULE PROCEDURE marching_init
        END INTERFACE
        
        
      CONTAINS
        
#include "marching_new.F90"
#include "marching_finalize.F90"
#include "marching_relax.F90"
#include "marching_marching.F90"
#include "marching_integrate.F90"
#include "marching_adjust_flow_v.F90"

      END MODULE Class_Marching
      
