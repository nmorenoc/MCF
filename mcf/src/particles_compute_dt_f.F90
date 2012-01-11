      SUBROUTINE particles_compute_dt_f(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_compute_dt_f
        !----------------------------------------------------
        !
        ! Purpose     : Adapt dt according to the new
        !               maximum accerlation of particles,
        !               compute dt.
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
        
        TYPE(Particles), INTENT(INOUT)  :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        
        !----------------------------------------------------
        ! Calculate dt according to maximum acceleration 
        ! constraints.
        !----------------------------------------------------
        
        IF( this%fa_max > 0.0_MK ) THEN
           
           this%dt_f = 0.25_MK * SQRT(this%h / this%fa_max)
           
        ELSE
           
           this%dt_f = -1.0_MK
           
        END IF
        
        RETURN
        
      END SUBROUTINE particles_compute_dt_f
      
      
