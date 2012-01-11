!------------------------------------------------------------
! All the public "set" subroutines of Class StateEquation,
! which return member variables of StateEquation.
!
! Reference   :
!
! Remark      :
!
! Revisions   :  V0.1 15.03.2010, original version.
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

      SUBROUTINE stateEquation_set_rho_ref(this,d_rho_ref,stat_info)
        
        TYPE(StateEquation),INTENT(OUT) :: this
        REAL(MK), INTENT(IN)            :: d_rho_ref
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%rho_ref = d_rho_ref
        
        RETURN
        
      END SUBROUTINE stateEquation_set_rho_ref
