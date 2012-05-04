!----------------------------------------------------------------
! This contains the routines which calculate the force and 
! return two accelerations between two Newtonian particles, 
! accelerations from i.e. conservative(pressure) force 
! dissipative(viscous) force and random(thermal noise) force.
! Currently, only Hu-Adams formulation has random force.
!----------------------------------------------------------------
#ifdef __PARTICLES_FORCE_SEPARATE
      SUBROUTINE rhs_force_ff_Newtonian(this,&
           xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&      
           mi,mj,w,gradw,fi,fj,fpi,fpj,fvi,fvj,fri,frj,auij,stat_info)
#else
     SUBROUTINE rhs_force_ff_Newtonian(this,&
           xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&      
           mi,mj,w,gradw,fi,fj,auij,stat_info)
#endif 
        !----------------------------------------------------
        ! Subroutine  :  Return pair-wise force per unit
        !                mass, i.e., acceleration between
        !                two Newtonian fluid particles.
        !                Currently, even we have colloids
        !                boundary particles, 
        !                fluid-boundary particles force
        !                calculation is treated same way.
        ! 
        !  Revision   : V0.2 10.08.2009
        !               Return two accleration of two particles
        !               instead of one, then it is able
        !               to handel particles with different
        !               masses.     
        !              
        !               V0.1 01.03.2009, original version
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
        !  fi        : force acting on i (per unit mass)
        !  fj        : force acting on j (per unit mass)
        !  fpi       : pressure force acting on i (per unit mass)
        !  fpj       : pressure force acting on j (per unit mass)
        !  fvi       : viscous force acting on i (per unit mass)
        !  fvj       : viscous force acting on j (per unit mass)
        !  fri       : random force acting on i (per unit mass)
        !  frj       : random force acting on j (per unit mass)
        !  auij      : potential energy between i,j
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
        REAL(MK), INTENT(IN)                    :: pi
        REAL(MK), INTENT(IN)                    :: pj
        REAL(MK), INTENT(IN)                    :: mi
        REAL(MK), INTENT(IN)                    :: mj
        REAL(MK), INTENT(IN)                    :: w
        REAL(MK), INTENT(IN)                    :: gradw
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: fi
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: fj
#ifdef __PARTICLES_FORCE_SEPARATE
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: fpi
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: fpj
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: fvi
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: fvj
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: fri
        REAL(MK), DIMENSION(:), INTENT(OUT), OPTIONAL     :: frj
#endif
        REAL(MK), INTENT(OUT), OPTIONAL         :: auij 
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
        
        
        IF(PRESENT(auij))THEN
           auij = 0.0_MK
        END IF
        
        
        !----------------------------------------------------
        ! Select different right hand side formulation.
        !----------------------------------------------------
        
        SELECT CASE (this%rhs_force_type)
           
        CASE (0) 
           
           PRINT *, "rhs_force_ff_Newtonian : ", &
                "0 not existing formulation !"
           stat_info = -1
           GOTO 9999
           
        CASE (1)
           
           CALL  rhs_force_ff_Newtonian_Morris(this,&
                xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&
                mi,mj,w,gradw,fi,fj,auij,stat_info_sub)
           
        CASE (2)
           
           CALL  rhs_force_ff_Newtonian_Espanol(this,&
                xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&
                mi,mj,w,gradw,fi,fj,auij,stat_info_sub)
           
        CASE (3)
           
           CALL  rhs_force_ff_Newtonian_HuAdams(this,&
                xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&
                mi,mj,w,gradw,fi,fj,auij,stat_info_sub)
           
        CASE (4)
           
#if __PARTICLES_FORCE_SEPARATE
           IF ( PRESENT(fpi) .AND. PRESENT (fpj) .AND. &
                PRESENT(fvi) .AND. PRESENT (fvj) .AND. &
                PRESENT(fri) .AND. PRESENT (frj) ) THEN
              
              CALL rhs_force_ff_Newtonian_HuAdams_angular(this,&
                   xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&
                   mi,mj,w,gradw,fi,fj,&
                   fpi,fpj,fvi,fvj,fri,frj,auij=auij,stat_info=stat_info_sub)
           ELSE
              
              CALL rhs_force_ff_Newtonian_HuAdams_angular(this,&
                   xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&
                   mi,mj,w,gradw,fi,fj,auij=auij,stat_info=stat_info_sub)
              
           END IF
#else
           CALL  rhs_force_ff_Newtonian_HuAdams_angular(this,&
                xi,xj,dij,vi,vj,rhoi,rhoj,pi,pj,&
                mi,mj,w,gradw,fi,fj,auij=auij,stat_info=stat_info_sub)
#endif     
           
        CASE DEFAULT
           
           PRINT *, "rhs_force_ff_Newtonian : ", & 
                "4-> Not existing formulation!"
           stat_info = -1
           GOTO 9999
           
        END SELECT ! rhs_force_type
        
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "rhs_force_ff_Newtonian :", &
                "rhs_force_ff_formulation has some problem!"
           stat_info = -1
           GOTO 9999
        END IF
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE rhs_force_ff_Newtonian
      
#include "rhs_force_ff_Newtonian_Morris.F90"    
#include "rhs_force_ff_Newtonian_Espanol.F90"
#include "rhs_force_ff_Newtonian_HuAdams.F90"
#include "rhs_force_ff_Newtonian_HuAdams_angular.F90"
