      SUBROUTINE colloid_noslip(this,xf,xc,vf,vc,sid_c,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_noslip
        !----------------------------------------------------
        !
        ! Purpose     : Return an artificial velocity for a 
        !               numerical particle inside a colloid
        !               object, in order to get no slip 
        !               condition on the surface of a colloid.
        !
        !               It will use accordingly, different
        !               ways of implementing no slip
        !               boundary condition.
        !
        !  Revision   : V0.1  01.03.2009, original version.
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
        ! Arguments :
        !
        ! Input 
        !
        ! this  : object of a colloid.
        ! xf    : position of a fluid particle.
        ! xc    : position of a colloid boundary particle.
        ! vf    : velocity of a fluid particle.
        ! sid_c : species ID of a colloid boundary particle.
        !
        ! Output
        !
        ! vc    : extrapolated velocity for the colloid
        !         boundary particle.
        ! stat_info : status of the routine.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xf
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xc
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vf
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: vc
        INTEGER, INTENT(IN)                     :: sid_c        
        INTEGER, INTENT(OUT)                    :: stat_info 
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub

        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        vc(:) = 0.0_MK
        
        !----------------------------------------------------
        ! Decide which no slip type to choose.
        !----------------------------------------------------
        
        SELECT CASE( this%noslip_type )
           
        CASE ( mcf_no_slip_frozen )
           
           CALL colloid_noslip_frozen(this,xc,vc,sid_c,stat_info_sub)
           
           IF( stat_info_sub /= 0 ) THEN
              PRINT *, "colloid_noslip : ", &
                   "Freezing boundary particle has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
        CASE ( mcf_no_slip_Morris )
           
           CALL colloid_noslip_Morris(this,xf,xc,vf,vc,sid_c,&
                stat_info_sub)
           
           IF( stat_info_sub /= 0 ) THEN
              PRINT *, "colloid_noslip : ", &
                   "Morris boundary particle has problem ! "
              stat_info = -1
              GOTO 9999
           END IF
           
        CASE (3)
           
           CALL colloid_noslip_Zhu(this,xf,xc,vf,vc,sid_c,&
                stat_info_sub)
           
           IF( stat_info_sub /= 0 ) THEN
              PRINT *, "colloid_noslip : ", &
                   "Zhu boundary particle has problem ! "
              stat_info = -1
              GOTO 9999
           END IF
           
        END SELECT
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_noslip
      
#include "colloid_noslip_frozen.F90"
#include "colloid_noslip_Morris.F90"
#include "colloid_noslip_Zhu.F90"

