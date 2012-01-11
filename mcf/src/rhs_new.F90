      SUBROUTINE rhs_init(this,d_ctrl,d_phys,d_random,stat_info)
        !----------------------------------------------------
        ! Subroutine  : rhs_init
        !----------------------------------------------------
        !
        ! Purpose     : Construtor of rhs Class.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 15.07.2009, original version.
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
        
        TYPE(Rhs), INTENT(OUT)                  :: this
        TYPE(Control), INTENT(IN), TARGET       :: d_ctrl
        TYPE(Physics), INTENT(IN), TARGET       :: d_phys
        TYPE(Random), INTENT(IN), TARGET        :: d_random
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------

        INTEGER                                 :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization.
        !----------------------------------------------------

        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(this%ctrl)
        NULLIFY(this%phys)
        NULLIFY(this%random)
        
        this%ctrl  => d_ctrl
        this%rhs_density_type =  &
             control_get_rhs_density_type(d_ctrl,stat_info_sub)
        this%Newtonian        = &
             control_get_Newtonian(d_ctrl,stat_info_sub)
        this%Brownian         = &
             control_get_Brownian(d_ctrl,stat_info_sub)
        this%rhs_force_type   =  &
             control_get_rhs_force_type(d_ctrl,stat_info_sub)
        this%phys   => d_phys
        this%num_dim = &
             physics_get_num_dim(d_phys,stat_info_sub)
        this%dt      = &
             physics_get_dt(d_phys,stat_info_sub)
        this%eta      = &
             physics_get_eta(d_phys,stat_info_sub)
        this%eta_coef = &
             physics_get_eta_coef(d_phys,stat_info_sub)
        this%kt      = &
             physics_get_kt(d_phys,stat_info_sub)
        
        this%random => d_random
        
        RETURN
        
      END SUBROUTINE rhs_init


      SUBROUTINE rhs_display_parameters(this,stat_info)

        TYPE(Rhs),INTENT(IN)          ::this
        INTEGER,INTENT(OUT)           :: stat_info
        
        stat_info = 0
        
        PRINT *, '------------------Start------------------'
        PRINT *, '     Rhs parameters'
        PRINT *, '-----------------------------------------'
        
        
        SELECT CASE(this%rhs_density_type) 
           
        CASE (1)
           
           PRINT *, "rhs_density_type : ", &
                "Mass density summation"
           
        CASE (2) 
           
           PRINT *, "rhs_density_type : ", &
                "Number density summation"
           
        CASE (3)
           
           PRINT *, "rhs_density_type :", &
                "Time evolution of mass density(not available yet)"
           stat_info = -1
           GOTO 9999
           
        END SELECT

        IF (this%Newtonian)  THEN
           
           PRINT *, "Newtonian        : ", &
                "Yes."
        ELSE
           
           PRINT *, "Newtonian        : ", &
                "No."
           
        END IF
      
        IF (this%Brownian)  THEN
           
           PRINT *, "Brownian         : ", &
                "Yes."
        ELSE
           
           PRINT *, "Brownian         : ", &
                "No."
           
        END IF
        
        SELECT CASE(this%rhs_force_type) 
           
        CASE (1)
           
           PRINT *, "rhs_force_type   : ", &
                "Morris et al., J. Comput. Phys. 1997"
           
        CASE (2) 
           
           PRINT *, "rhs_force_type   : ", &
                "Espanol and Revenga, Phys. Rev. E 2003"
           
        CASE (3) 
           
           PRINT *, "rhs_force_type   : ", &
                "Hu and Adams, J. Comput. Phys. 2006"
           
        CASE (4) 
           
           PRINT *, "rhs_force_type   : ", &
                "Hu and Adams, Phys. Fluids. 2006"
       
        END SELECT

        
        PRINT *, '-------------------End-------------------'
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE rhs_display_parameters
