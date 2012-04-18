      SUBROUTINE marching_integrate(this,step,time,dt,stat_info)
        !----------------------------------------------------
        ! SUBROUTINE  : marching_integrate
        !
        ! Purpose     : Integrator control routine/
        !               interface called by outside.
        !               The desired intergration scheme
        !               will be used here correctly.
        !
        ! Revision    : V0.1 01.04 2009, original version.
        !
        !----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        ! 
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
      
         
        !----------------------------------------------------
        ! Arguments
        ! this      : Marching object
        ! step      : current step.
        ! time      : current time.
        ! dt        : time step.
        ! stat_info : return status.
        !----------------------------------------------------
        
        TYPE(Marching), INTENT(INOUT)   :: this
        INTEGER, INTENT(IN)             :: step
        REAL(MK), INTENT(IN)            :: time
        REAL(MK), INTENT(IN)            :: dt
        INTEGER, INTENT(OUT)	        :: stat_info
        
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: integrate_type
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        integrate_type = &
             control_get_integrate_type(this%ctrl,stat_info_sub)
        
        SELECT CASE(integrate_type) 
           
        CASE (1)
           
           CALL marching_integrate_Euler(this,step,time,dt,stat_info_sub)
           
        CASE (2)
           
           CALL marching_integrate_VV(this,step,time,dt,stat_info_sub)
           
        CASE DEFAULT
           
           PRINT *, "marching_integrate: ", &
                "integrator not available! "
           stat_info = -1
           GOTO 9999
           
        END SELECT
        
        IF( stat_info_sub /=0 ) THEN
           PRINT *, "marching_integrate: ", &
                "Integrating failed !" 
           stat_info = -1
           GOTO 9999
        END IF
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE marching_integrate
      
#include "marching_integrate_Euler.F90"
#include "marching_integrate_VV.F90"
