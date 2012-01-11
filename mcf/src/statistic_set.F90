!------------------------------------------------------------
! All the public "set" subroutines of Class statistic,
! which return member variables of statistic.
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

      SUBROUTINE statistic_set_rho_min(this,d_rho_min,stat_info)
        
        TYPE(Statistic), INTENT(OUT)    :: this
        REAL(MK), INTENT(IN)            :: d_rho_min
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        
        this%rho_min = d_rho_min
        
        RETURN
        
      END SUBROUTINE statistic_set_rho_min
      

      SUBROUTINE statistic_set_rho_max(this,d_rho_max,stat_info)
        
        TYPE(Statistic), INTENT(OUT)    :: this
        REAL(MK), INTENT(IN)            :: d_rho_max
        INTEGER, INTENT(OUT)            :: stat_info 
        
        stat_info = 0
        
        this%rho_max = d_rho_max
        
        RETURN        
        
      END SUBROUTINE statistic_set_rho_max
      
