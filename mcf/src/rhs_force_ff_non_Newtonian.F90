!----------------------------------------------------------------
! This contains the routines which calculate the force between
! two Non-Newtonian particles, i.e. conservative(pressure) force 
! dissipative(viscous) force and random(thermal noise) force.
! Currently, only HuAdams formulation has random force.
!----------------------------------------------------------------
      SUBROUTINE rhs_force_ff_non_Newtonian(this,&
           xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&      
           mi,mj,w,gradw,fi,fj,stat_info)
        
        !----------------------------------------------------
        ! Subroutine  :  Return pair-wise force between
        !                two Non-Newtonian fluid particles.
        !                Currently, even we have colloids,
        !                fluid-colloid force is treated same
        !                way.
        !                Now we have only Hu_Admas formulation
        !                for Non-Newtonian.
        ! 
        ! Revision    :  V0.1 09.03.2011
        !                The non Newtonian force for the 
        !                Espanol model is implimented.(Adolfo)
        !
        !                V0.1 25.08.2010
        !                The dt variable was not calculated.
        !                Now it is. (Adolfo)
        !                
        !                V0.1 01.03.2009
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
        
        !----------------------------------------------------
        !  Arguments
        !  xi        : particle i 's position
        !  xj        : particle j 's position
        !  dij       : distance between i,j
        !  vi        : i's velocity
        !  vj        : j's velocity
        !  rhoi      : i's density
        !  rhoj      : j's density
        !  pi        : i's pressure
        !  pj        : j's pressure
        !  mi        : i's mass
        !  mj        : j's mass
        !  w         : kernel value
        !  gradw     : kernel gradient value
        !  fij       : force between i,j
        !  stat_info : return flag of status
        !----------------------------------------------------
        
        TYPE(Rhs), INTENT(INOUT)                :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xi
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xj
        REAL(MK), INTENT(IN)                    :: dij
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vi
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vj
        REAL(MK), INTENT(IN)                    :: rhoi
        REAL(MK), INTENT(IN)                    :: rhoj
        REAL(MK), DIMENSION(:,:),INTENT(IN)     :: pi
        REAL(MK), DIMENSION(:,:),INTENT(IN)     :: pj
        REAL(MK), INTENT(IN)                    :: mi
        REAL(MK), INTENT(IN)                    :: mj
        REAL(MK), INTENT(IN)                    :: w
        REAL(MK), INTENT(IN)                    :: gradw
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: fi
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: fj
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
        

        !----------------------------------------------------
        ! Select different right hand side formulation.
        !----------------------------------------------------
        
        SELECT CASE (this%rhs_force_type)
           
        CASE (2)
           
           CALL  rhs_force_ff_non_Newtonian_Espanol(this,&
                xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&
                mi,mj,w,gradw,fi,fj,stat_info_sub)

        CASE (3)
           
           CALL  rhs_force_ff_non_Newtonian_HuAdams(this,&
                xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&
                mi,mj,w,gradw,fi,fj,stat_info_sub)
           
        case (4)

           CALL  rhs_force_ff_non_Newtonian_HuAdams_angular(this,&
                xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&
                mi,mj,w,gradw,fi,fj,stat_info_sub)
           
        CASE DEFAULT
           
           PRINT *, "rhs_force_ff_non_Newtonian : ", & 
                "Not existing formulation !"
           stat_info = -1
           GOTO 9999
           
        END SELECT ! rhs_force_type
        
        IF(stat_info_sub /= 0) THEN
           PRINT *, "rhs_force_ff_non_Newtonian :", &
                "formulation has some problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE rhs_force_ff_non_Newtonian
      

#include "rhs_force_ff_non_Newtonian_Espanol.F90"
#include "rhs_force_ff_non_Newtonian_HuAdams.F90"
#include "rhs_force_ff_non_Newtonian_HuAdams_angular.F90"
      
