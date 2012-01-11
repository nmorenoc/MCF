      SUBROUTINE stateEquation_compute_pressure_scalar(this,rho,p,stat_info)
        !----------------------------------------------------
      	! Subroutine : StateEquation_compute_pressure_scalar
	!----------------------------------------------------
      	!
      	! Purpose    :
	!> \brief      Computing pressure scalar variable
        !>             using state equation.
      	!>	   	
        !
        ! References :
     	!
      	! Revisions  : V0.1 03.03.2009, original version.
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
        
        TYPE(StateEquation), INTENT(IN)         :: this
        REAL(MK), INTENT(IN)                    :: rho
        REAL(MK), INTENT(OUT)                   :: p
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: stat_info_sub
        
        
        stat_info     = 0
        stat_info_sub = 0

        SELECT CASE (this%stateEquation_type)
           
        CASE (1)
           
           CALL stateEquation_compute_pressure_scalar_Morris(this,rho,p,stat_info_sub)
           
        CASE (2)
           
           CALL stateEquation_compute_pressure_scalar_Batchelor(this,rho,p,stat_info_sub)   
           
        END SELECT
        
        RETURN        
        
      END SUBROUTINE stateEquation_compute_pressure_scalar

      
      SUBROUTINE stateEquation_compute_pressure_scalar_Morris(this,rho,p,stat_info)

        TYPE(StateEquation), INTENT(IN)         :: this
        REAL(MK), INTENT(IN)                    :: rho
        REAL(MK), INTENT(OUT)                   :: p
        INTEGER, INTENT(OUT)                    :: stat_info


        stat_info = 0

        p = (this%c)**2 * (rho-this%rho_ref)

        RETURN

      END SUBROUTINE stateEquation_compute_pressure_scalar_Morris


      SUBROUTINE stateEquation_compute_pressure_scalar_Batchelor(this,rho,p,stat_info)

        TYPE(StateEquation), INTENT(IN)         :: this
        REAL(MK), INTENT(IN)                    :: rho
        REAL(MK), INTENT(OUT)                   :: p
        INTEGER, INTENT(OUT)                    :: stat_info


        stat_info = 0

        p = this%p0 * ( (rho/this%rho_ref)**this%gamma -1 )

        RETURN
        
      END SUBROUTINE stateEquation_compute_pressure_scalar_Batchelor


      SUBROUTINE stateEquation_compute_pressure_vector(this,rho,p,num_part,stat_info)
         !----------------------------------------------------
      	! Subroutine : StateEquation_compute_pressure_scalar
	!----------------------------------------------------
      	!
      	! Purpose    :
	!> \brief      Computing pressure vector variables
        !>             using state equation.
      	!>	   	
        !
        ! References :
     	!
      	! Revisions  : V0.1 03.03.2009, original version.
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
        
        TYPE(StateEquation), INTENT(IN)         :: this
        REAL(MK), DIMENSION(:),INTENT(IN)       :: rho
        REAL(MK), DIMENSION(:), POINTER         :: p
        INTEGER, INTENT(IN)                     :: num_part
        INTEGER, INTENT(OUT)                    :: stat_info

        INTEGER                                 :: stat_info_sub

        stat_info     = 0
        stat_info_sub = 0

        IF(ASSOCIATED(p)) THEN
           DEALLOCATE(p)
        END IF
        
        ALLOCATE(p(num_part))
        
        SELECT CASE (this%stateEquation_type)
           
        CASE (1)
           
           CALL stateEquation_compute_pressure_vector_Morris(this, &
                rho,p,num_part,stat_info_sub)
           
        CASE (2)
           
           CALL stateEquation_compute_pressure_vector_Batchelor(this, &
                rho,p,num_part,stat_info_sub)
           
        END SELECT
      
        RETURN        
        
      END SUBROUTINE stateEquation_compute_pressure_vector

      
      SUBROUTINE stateEquation_compute_pressure_vector_Morris(this,rho,p,num_part,stat_info)
        
        TYPE(StateEquation), INTENT(IN)         :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: rho
        REAL(MK), DIMENSION(:), POINTER         :: p
        INTEGER, INTENT(IN)                     :: num_part
        INTEGER, INTENT(OUT)                    :: stat_info


        stat_info = 0

        p(1:num_part) = (this%c)**2 * ( rho(1:num_part)**this%gamma-this%rho_ref)

        RETURN

      END SUBROUTINE stateEquation_compute_pressure_vector_Morris

      
      SUBROUTINE stateEquation_compute_pressure_vector_Batchelor(this,rho,p,num_part,stat_info)

        TYPE(StateEquation), INTENT(IN)         :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: rho
        REAL(MK), DIMENSION(:), POINTER         :: p
        INTEGER, INTENT(IN)                     :: num_part
        INTEGER, INTENT(OUT)                    :: stat_info        
        
        
        stat_info = 0

        p(1:num_part) = this%p0 * ( (rho(1:num_part)/this%rho_ref)**this%gamma -1 )

        RETURN
        
      END SUBROUTINE stateEquation_compute_pressure_vector_Batchelor
