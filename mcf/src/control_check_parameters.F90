      LOGICAL FUNCTION control_check_parameters(this,stat_info)
        !----------------------------------------------------
        !  Subroutine   :  control_check_parameters
        !----------------------------------------------------
        !
        !  Purpose      :  Check if control parameters are
        !                  resonable.
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : V0.1 15.07.2009, original version.
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
        TYPE(Control), INTENT(INOUT)    :: this
        INTEGER, INTENT(OUT)            :: stat_info

        
        stat_info = 0
        
        control_check_parameters = .TRUE.

        IF ( (this%Brownian .EQV. .TRUE.) .AND. &
             ( this%symmetry .NEQV. .TRUE.) ) THEN
           
           control_check_parameters = .FALSE.
           PRINT *, "control_check_parameters : ", &
                "For Brownian solvent, do use symmetry interaction and &
                symmetry inter-processor communication !"
           GOTO 9999
           
        END IF

9999    CONTINUE  
        
        RETURN        
        
      END FUNCTION control_check_parameters
      
