      MODULE Class_Tool
        !----------------------------------------------------
      	!       Class_Tool
      	!----------------------------------------------------
        !
      	!  Purpos     : auxiliary functions for global usage.
      	!
      	!  Remarks    :
      	!
      	!  References :
      	!
      	!  Version    : V 0.1. 03.03.2009, original version.
      	!----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
   
        USE mcf_header
        
        IMPLICIT NONE
        SAVE
        
        TYPE Tool
           INTEGER      :: flag
        END TYPE Tool
        
        
        INTERFACE tool_new
           MODULE PROCEDURE tool_init           
        END INTERFACE
        
      CONTAINS
        
#include "tool_new.F90"
#include "tool_uppercase.F90"
#include "tool_cross_product.F90"
#include "tool_rotation_matrix.F90"
#include "tool_rotation_vector.F90"
#include "tool_solve_linear_equations.F90"
#include "tool_L1_norm.F90"
#include "tool_L2_norm.F90"
#include "tool_Linfty_norm.F90"

      END MODULE Class_Tool
