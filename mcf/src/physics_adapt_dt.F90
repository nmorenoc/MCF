      SUBROUTINE physics_adapt_dt(this,dt_f,stat_info)
        !----------------------------------------------------
        ! Subroutine  : physics_adapt_dt
        !----------------------------------------------------
        !
        ! Purpose     : Adapt dt according to the new
        !               maximum accerlation of particles.
        !      
        ! Reference   : Morris et al. JCP 1997.
        !
        ! Remark      :
     
        !
        ! Revisions   : V0.1 13.10.20101, original version.
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
        !----------------------------------------------------
        
        TYPE(Physics), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: dt_f
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
       
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info  = 0
        
        this%dt_f  = dt_f
        
        
        !-------------------------------------------------
        ! Suppose to be positive infinity here.
        !-------------------------------------------------
        
        this%dt = 1.0e+6_MK
        
        IF ( this%dt_c > 0.0_MK .AND. &
             this%dt_c < this%dt )  THEN
           
           this%dt = this%dt_c
           
        END IF
        
        IF ( this%dt_nu > 0.0_MK .AND. &
             this%dt_nu < this%dt ) THEN
           
           this%dt = this%dt_nu
           
        END IF
        
        IF ( this%dt_f > 0.0_MK .AND. &
             this%dt_f < this%dt )  THEN
           
           this%dt = this%dt_f
           
        END IF
        
        
        RETURN
        
      END SUBROUTINE physics_adapt_dt
      
      
