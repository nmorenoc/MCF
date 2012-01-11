      MODULE Class_Debug
        !----------------------------------------------------
      	!  Class      :	Debug
	!----------------------------------------------------
      	!
      	!  Purpose    :
	!> \brief       Variables and corresponding operations
        !>              for Debug quantities.
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
        
        USE ppm_module_time
        USE mcf_header
        
        IMPLICIT NONE
        SAVE
        
        TYPE Debug
           PRIVATE
           INTEGER                      :: flag
           CHARACTER(LEN=MAX_CHAR)      :: time_file
           INTEGER                      :: time_file_unit
        END TYPE Debug
        
        TYPE(Debug)                     :: global_debug

        INTERFACE debug_new
           MODULE PROCEDURE debug_init
        END INTERFACE

        INTERFACE debug_substart
           MODULE PROCEDURE debug_substart_p
           MODULE PROCEDURE debug_substart_s
        END INTERFACE

        INTERFACE debug_substop
           MODULE PROCEDURE debug_substop_p
           MODULE PROCEDURE debug_substop_s
        END INTERFACE

        INTERFACE debug_print_msg
           MODULE PROCEDURE  debug_print_msg_a
           MODULE PROCEDURE  debug_print_msg_i
           MODULE PROCEDURE  debug_print_msg_f
           MODULE PROCEDURE  debug_print_msg_aa
           MODULE PROCEDURE  debug_print_msg_ai
           MODULE PROCEDURE  debug_print_msg_af
        END INTERFACE

        INTERFACE debug_write_output
           MODULE PROCEDURE debug_write_output_1d_f
           MODULE PROCEDURE debug_write_output_1d_i
           MODULE PROCEDURE debug_write_output_2d_f
           MODULE PROCEDURE debug_write_output_2d_i
        END INTERFACE

      CONTAINS
            
#include "debug_new.F90"
#include "debug_finalize.F90"
#include "debug_substart.F90"
#include "debug_substop.F90"
#include "debug_print_msg.F90"
#include "debug_write_output.F90"
#include "debug_validate_motion.F90"
#include "debug_get.F90"
#include "debug_set.F90"
#include "debug_open.F90"
#include "debug_close.F90"
#include "debug_write_time.F90"

      END MODULE Class_Debug
