!--------------------------------------------------
!  Subroutine   :  rhs_set_*
!--------------------------------------------------
!
!  Purpose      : Set routines of Class RHS.
!
!  Reference    :
!
!  Remark       :
!
!  Revisions    : V0.1 01.03.2009, original version.
!
!--------------------------------------------------
! Author       : Xin Bian
! Contact      : xin.bian@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!--------------------------------------------------      

      SUBROUTINE rhs_set_Brownian(this,d_Brownian,stat_info)

        TYPE(Rhs), INTENT(INOUT)        :: this
        LOGICAL, INTENT(IN)             :: d_Brownian
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0
        
        this%Brownian = d_Brownian     
        
        RETURN
        
      END SUBROUTINE rhs_set_Brownian
      
      
      SUBROUTINE rhs_set_rhs_density_type(this,d_rhs_density_type,stat_info)
        
        TYPE(Rhs), INTENT(INOUT)        :: this
        INTEGER, INTENT(IN)             :: d_rhs_density_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%rhs_density_type = d_rhs_density_type        
        
        RETURN
        
      END SUBROUTINE rhs_set_rhs_density_type


      SUBROUTINE rhs_set_dt(this,d_dt,stat_info)
        
        TYPE(Rhs), INTENT(INOUT)        :: this
        REAL(MK), INTENT(IN)            :: d_dt
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0
        
        this%dt = d_dt
        
        RETURN
        
      END SUBROUTINE rhs_set_dt


      SUBROUTINE rhs_set_kt(this,d_kt,stat_info)
        
        TYPE(Rhs), INTENT(INOUT)        :: this
        REAL(MK), INTENT(IN)            :: d_kt
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0
        
        this%kt = d_kt
        
        RETURN
        
      END SUBROUTINE rhs_set_kt
