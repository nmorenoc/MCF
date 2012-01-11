      MODULE Class_Tool
      
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
        
      END MODULE Class_Tool
