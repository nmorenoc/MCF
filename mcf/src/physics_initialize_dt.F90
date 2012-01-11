      SUBROUTINE physics_initialize_dt(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : physics_initialize_dt
        !----------------------------------------------------
        !
        ! Purpose     : Give dt an initial value according
        !               to CFL, viscous constrains and
        !               maximum body force.
        !      
        ! Reference   : Morris et al. JCP 1997.
        
        !
        ! Remark      : If dt is given from input as
        !               positive value, it will not be
        !               overritten. However, dt_nu etc
        !               will be caculated anyway as
        !               references.
     
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
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
       
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: relax_run
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        
        relax_run    = &
             control_get_relax_run(this%ctrl,stat_info_sub)
        
        !----------------------------------------------------
        ! Set basic initial dts
        !----------------------------------------------------
        
        this%dt_c  = -1.0_MK
        this%dt_nu = -1.0_MK
        this%dt_f  = -1.0_MK
        this%dt_c_relax = -1.0_MK
        
        
        !----------------------------------------------------
        ! CFL condition
        !----------------------------------------------------
        
        IF ( this%c > 0 ) THEN
           
           this%dt_c = 0.25_MK * this%h / this%c 
           
        END IF
        
        IF ( relax_run .AND. this%c_relax > 0 ) THEN
           
           this%dt_c_relax = 0.25_MK * this%h / this%c_relax 
           
        END IF
        
        !----------------------------------------------------
        ! Constraints due to viscous diffusion
        !----------------------------------------------------
        
        IF( this%eta > 0.0_MK ) THEN
           
           this%dt_nu = 0.125_MK * &
                (this%h**2) * this%rho / this%eta
           
        END IF
        
        !----------------------------------------------------
        ! Constraints according to acceleration
        !----------------------------------------------------
        
        IF( this%fa_max > 0.0_MK ) THEN
           
           this%dt_f = 0.25_MK * SQRT(this%h / this%fa_max)
           
        END IF
        
        !----------------------------------------------------
        ! If only non-positive value is given from input,
        ! take the minimum of dt_c, dt_nu, dt_f.
        !----------------------------------------------------
        
        IF ( this%dt <= 0.0_MK ) THEN
           
           !-------------------------------------------------
           ! Suppose to be infinity here.
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
           
        END IF ! this%dt < = 0
        
        
        !----------------------------------------------------
        ! For relax run, if only non-positive value is 
        ! given from input, take the minium of
        ! dt_c_relax, dt_nu, dt_f.
        !----------------------------------------------------
        
        IF ( relax_run .AND. this%dt_relax <= 0.0_MK ) THEN
           
           !-------------------------------------------------
           ! Suppose to be infinity here.
           !-------------------------------------------------
           
           this%dt_relax = 1.0e+6_MK
           
           IF ( this%dt_c_relax > 0.0_MK .AND. &
                this%dt_c_relax < this%dt_relax )  THEN
              
              this%dt_relax = this%dt_c_relax
              
           END IF
           
           IF ( this%dt_nu > 0.0_MK .AND. &
                this%dt_nu < this%dt_relax ) THEN
              
              this%dt_relax = this%dt_nu
              
           END IF
           
           IF ( this%dt_f > 0.0_MK .AND. &
                this%dt_f < this%dt_relax )  THEN
              
                 this%dt_relax = this%dt_f
                 
           END IF
           
        END IF ! this%dt_relax < = 0
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE physics_initialize_dt
      
      
