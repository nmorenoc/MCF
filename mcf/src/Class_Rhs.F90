      MODULE Class_Rhs
        
        USE mcf_header
        USE Class_Debug
        
        USE Class_Physics
        USE Class_Random
        
        IMPLICIT NONE
        SAVE
        
        TYPE Rhs
           PRIVATE
           TYPE(Control), POINTER       :: ctrl
           INTEGER                      :: rhs_density_type
           LOGICAL                      :: Newtonian
           LOGICAL                      :: Brownian
           INTEGER                      :: rhs_force_type
           TYPE(Physics), POINTER       :: phys
           INTEGER                      :: num_dim
           REAL(MK)                     :: dt
           REAL(MK)                     :: eta
           REAL(MK)                     :: eta_coef
           REAL(MK)                     :: kt
           TYPE(Random), POINTER        :: random
        END TYPE Rhs
        
        INTERFACE rhs_new
           MODULE PROCEDURE rhs_init
        END INTERFACE
        
        INTERFACE rhs_force_ff
           MODULE PROCEDURE rhs_force_ff_Newtonian
           MODULE PROCEDURE rhs_force_ff_non_Newtonian
        END INTERFACE
        
        INTERFACE rhs_force_cc
           MODULE PROCEDURE rhs_force_cc_Newtonian
        END INTERFACE
        
      CONTAINS
        
#include "rhs_new.F90"
#include "rhs_finalize.F90"
#include "rhs_set.F90"
#include "rhs_get.F90"
#include "rhs_density.F90"
#include "rhs_force_ff_Newtonian.F90"
#include "rhs_force_ff_non_Newtonian.F90"
#include "rhs_force_cc_Newtonian.F90"

      END MODULE Class_Rhs
      
