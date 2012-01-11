      MODULE Class_IO 
        !----------------------------------------------------
      	! Class	      :	IO
	!----------------------------------------------------
      	!
      	!  Purpose    :
	!> \brief       Variables and corresponding operations
        !>              for IO infomation.
      	!>
        !>              The variable memebers are private
        !
	! Remarks     :
      	!
      	! References  :
     	!
      	! Revisions   : V0.2 08.11. 2010, introduce output
        !               frequency also by time(only step before).
        !               It is very useful especially for
        !               adaptive dt.
        !
        !               V0.1 03.03.2009
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
        USE Class_Tool
        
        USE Class_Control
        USE Class_Physics
        USE Class_Particles
        USE Class_Statistic
        
        IMPLICIT NONE
        SAVE
        
        TYPE IO
           PRIVATE
           
           !-------------------------------------------------
           ! ctrl : control object.
           ! phys : physics object.
           !-------------------------------------------------
           
           
           TYPE(Control), POINTER       :: ctrl
           TYPE(Physics), POINTER       :: phys
           
           !-------------------------------------------------
           ! write_output : Indicating if output is written.
           !                0 : no
           !                1 : according to step frequency
           !                2 : according to time frequency
           !
           !-------------------------------------------------
           
           INTEGER                      :: write_output
           INTEGER                      :: step_start
           REAL(MK)                     :: time_start
           
           !-------------------------------------------------
           ! *_file : name of the file.
           ! *_unit : UNIT for opening the file.
           ! *_fmt  : format of the file.
           ! *_freq : for writting, the frequency.
           !-------------------------------------------------
           
           CHARACTER(LEN=MAX_CHAR)      :: ctrl_file
           INTEGER                      :: ctrl_unit

           CHARACTER(LEN=MAX_CHAR)      :: physics_config_file
           INTEGER                      :: physics_config_unit
           
           CHARACTER(LEN=MAX_CHAR)      :: io_config_file
           INTEGER                      :: io_config_unit
           
           CHARACTER(LEN=MAX_CHAR)      :: read_particles_file
           INTEGER                      :: read_particles_unit
           CHARACTER(LEN=MAX_CHAR)      :: read_particles_fmt
           
           CHARACTER(LEN=MAX_CHAR)      :: read_conformation_file
           INTEGER                      :: read_conformation_unit
           CHARACTER(LEN=MAX_CHAR)      :: read_conformation_fmt
           
           CHARACTER(LEN=MAX_CHAR)      :: output_particles_relax_file
           INTEGER                      :: output_particles_relax_unit
           CHARACTER(LEN=MAX_CHAR)      :: output_particles_relax_fmt
           INTEGER                      :: output_particles_relax_freq_step

           CHARACTER(LEN=MAX_CHAR)      :: output_particles_file
           INTEGER                      :: output_particles_unit
           CHARACTER(LEN=MAX_CHAR)      :: output_particles_fmt
           INTEGER                      :: output_particles_freq_step
           REAL(MK)                     :: output_particles_freq_time
           INTEGER                      :: output_particles_freq_time_num

           CHARACTER(LEN=MAX_CHAR)      :: output_conformation_file
           INTEGER                      :: output_conformation_unit
           CHARACTER(LEN=MAX_CHAR)      :: output_conformation_fmt
           INTEGER                      :: output_conformation_freq_step
           REAL(MK)                     :: output_conformation_freq_time
           INTEGER                      :: output_conformation_freq_time_num

           CHARACTER(LEN=MAX_CHAR)      :: colloid_file
           INTEGER                      :: colloid_unit
           CHARACTER(LEN=MAX_CHAR)      :: colloid_fmt
           INTEGER                      :: colloid_freq_step
           REAL(MK)                     :: colloid_freq_time
           INTEGER                      :: colloid_freq_time_num

           CHARACTER(LEN=MAX_CHAR)      :: statistic_relax_file
           INTEGER                      :: statistic_relax_unit
           CHARACTER(LEN=MAX_CHAR)      :: statistic_relax_fmt
           INTEGER                      :: statistic_relax_freq_step
           
           CHARACTER(LEN=MAX_CHAR)      :: statistic_file
           INTEGER                      :: statistic_unit
           CHARACTER(LEN=MAX_CHAR)      :: statistic_fmt
           INTEGER                      :: statistic_freq_step
           REAL(MK)                     :: statistic_freq_time
           INTEGER                      :: statistic_freq_time_num

           CHARACTER(LEN=MAX_CHAR)      :: boundary_file
           INTEGER                      :: boundary_unit
           CHARACTER(LEN=MAX_CHAR)      :: boundary_fmt
           INTEGER                      :: boundary_freq_step
           REAL(MK)                     :: boundary_freq_time
           INTEGER                      :: boundary_freq_time_num
        
           CHARACTER(LEN=MAX_CHAR)      :: restart_particles_relax_file
           INTEGER                      :: restart_particles_relax_unit
           CHARACTER(LEN=MAX_CHAR)      :: restart_particles_relax_fmt

           CHARACTER(LEN=MAX_CHAR)      :: restart_physics_file
           INTEGER                      :: restart_physics_unit
           CHARACTER(LEN=MAX_CHAR)      :: restart_physics_fmt
           
           CHARACTER(LEN=MAX_CHAR)      :: restart_particles_file
           INTEGER                      :: restart_particles_unit
           CHARACTER(LEN=MAX_CHAR)      :: restart_particles_fmt
           
           CHARACTER(LEN=MAX_CHAR)      :: restart_conformation_file
           INTEGER                      :: restart_conformation_unit
           CHARACTER(LEN=MAX_CHAR)      :: restart_conformation_fmt
           
           INTEGER                      :: restart_freq_step
           REAL(MK)                     :: restart_freq_time
           REAL(MK)                     :: restart_freq_time_wall
           INTEGER                      :: restart_freq_time_num
           
           !-------------------------------------------------
           ! Flags indicating whether we should
           ! write corresponding file.
           !-------------------------------------------------
           
           INTEGER                      :: write_restart
           LOGICAL                      :: write_particles
           LOGICAL                      :: write_conformation
           LOGICAL                      :: write_colloid           
           LOGICAL                      :: write_statistic
           LOGICAL                      :: write_boundary
           LOGICAL                      :: write_restart_physics
           LOGICAL                      :: write_restart_particles
           LOGICAL                      :: write_restart_conformation
           
           !-------------------------------------------------
           ! Object of an auxiliary tool.
           !-------------------------------------------------
           
           TYPE(Tool)                   :: io_tool
           
        END TYPE IO
        
        
        INTERFACE io_new
           MODULE PROCEDURE io_init_default
           MODULE PROCEDURE io_init
        END INTERFACE
        
      CONTAINS
        
#include "io_new.F90"
#include "io_finalize.F90"
#include "io_check_parameters.F90"
#include "io_adjust_parameters.F90"
#include "io_get.F90"
#include "io_open.F90"
#include "io_read.F90"
#include "io_write.F90"
#include "io_write_condition.F90"
#include "io_close.F90"
        
      END MODULE Class_IO
      
