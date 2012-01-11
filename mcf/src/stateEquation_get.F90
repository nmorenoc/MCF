!------------------------------------------------------------
! All the public "get" subroutines of Class StateEquation,
! which return member variables of StateEquation.
!
! Reference   :
!
! Remark      :
!
! Revisions   :  V0.1 03.03.2009, original version.
!
!------------------------------------------------------------
! Author       : Xin Bian
! Contact      : xin.bian@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!------------------------------------------------------------

      REAL(MK) FUNCTION stateEquation_get_c(this,stat_info)

        TYPE(StateEquation), INTENT(IN) :: this
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0
        stateEquation_get_c = this%c

        RETURN
        
      END FUNCTION stateEquation_get_c

      
      REAL(MK) FUNCTION stateEquation_get_p0(this,stat_info)

        TYPE(StateEquation), INTENT(IN) :: this
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0
        stateEquation_get_p0 = this%p0

        RETURN
        
      END FUNCTION stateEquation_get_p0

      
      REAL(MK) FUNCTION stateEquation_get_rho_ref(this,stat_info)
        
        TYPE(StateEquation), INTENT(IN) :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        stateEquation_get_rho_ref = this%rho_ref

        RETURN
        
      END FUNCTION stateEquation_get_rho_ref

      
      REAL(MK) FUNCTION stateEquation_get_gamma(this,stat_info)
        
        TYPE(StateEquation), INTENT(IN) :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        stateEquation_get_gamma = this%gamma

        RETURN
        
      END FUNCTION stateEquation_get_gamma
