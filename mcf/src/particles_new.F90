      SUBROUTINE particles_init(this,d_ctrl,d_phys,&
           d_rhs,d_stateEquation,d_kern,d_tech,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_init
        !----------------------------------------------------
        !
        ! Purpose     :   Constructor of Particles object,
        !                  for initialization of its variables.
        !
        ! Input       :
        !
        ! Input/output: 
        !
        ! Output      : stat_info     (I) return status
        !
        ! Routines    :
        !
        ! Remarks     :
        !
        ! References  :
        !
        ! Revisions   : V0.1 01.03.2009, original version.
        !----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.     
        !----------------------------------------------------
        
        TYPE(Particles),INTENT(OUT)             :: this
        TYPE(Control), INTENT(IN),TARGET        :: d_ctrl  
        TYPE(Physics), INTENT(IN),TARGET        :: d_phys
        TYPE(Rhs), INTENT(IN), TARGET           :: d_rhs
        TYPE(StateEquation),INTENT(IN),TARGET   :: d_stateEquation
        TYPE(Kernel),INTENT(IN),TARGET          :: d_kern
        TYPE(Technique),INTENT(IN),TARGET       :: d_tech
        INTEGER,INTENT(OUT)                     :: stat_info

        !----------------------------------------------------
        !Local variables
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num_colloid 
        TYPE(Colloid),POINTER                   :: colloids
        
        stat_info     = 0
        NULLIFY(colloids)
        
        this%ctrl    => d_ctrl
        this%phys    => d_phys
        this%rhs     => d_rhs
        this%stateEquation => d_stateEquation
        this%tech    => d_tech
        this%kern    => d_kern
        
        this%pp_interact_cc = &
             control_get_pp_interact_cc(this%ctrl,stat_info_sub)
        this%pp_interact_cw = &
             control_get_pp_interact_cw(this%ctrl,stat_info_sub)
        
        this%num_dim = &
             physics_get_num_dim(this%phys,stat_info_sub)
        this%h       = &
             physics_get_h(this%phys,stat_info_sub)
        NULLIFY(this%x)
        NULLIFY(this%v)
        NULLIFY(this%rho)
        NULLIFY(this%rho_norm)
        this%rho_min  = &
             physics_get_rho(this%phys,stat_info_sub)
        this%rho_max  = &
             physics_get_rho(this%phys,stat_info_sub)
        NULLIFY(this%m)
        NULLIFY(this%p)
        NULLIFY(this%id)
        this%num_id   = 2
        this%pid_idx  = 1
        this%sid_idx  = 2         
        NULLIFY(this%f)
        this%fa_min = 0.0_MK
        this%fa_max = 0.0_MK

        NULLIFY(this%u)
        NULLIFY(this%au)

        NULLIFY(this%vgt)
        
        NULLIFY(this%evgt)
        NULLIFY(this%eval)
        NULLIFY(this%aeval)
        NULLIFY(this%evec)
        NULLIFY(this%aevec)

        NULLIFY(this%ct)
        NULLIFY(this%act)
        
        NULLIFY(this%pt)

        this%num_part_real   = 0
        this%num_part_all    = 0
        this%num_part_ghost  = 0
        
        this%num_part_fluid      = 0
        this%num_part_sym        = 0
        this%num_part_wall_sym   = 0
        this%num_part_wall_solid = 0
        this%num_part_wall_solid_real  = 0
        this%num_part_wall_solid_ghost = 0
        this%num_part_le         = 0
        
        NULLIFY(this%part_sym_list)
        NULLIFY(this%part_wall_sym_list)
        NULLIFY(this%part_wall_solid_real_list)
        NULLIFY(this%part_wall_solid_ghost_list)
        NULLIFY(this%part_le_list)
        
        this%num_part_colloid = 0
        NULLIFY(this%part_colloid_list)
        
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF ( num_colloid > 0 ) THEN
           
           CALL physics_get_colloid(this%phys,colloids,stat_info_sub)
           CALL colloid_set_tech(colloids,this%tech,stat_info_sub)
           
        END IF
        
        CALL tool_new(this%tool,stat_info_sub)
        
        NULLIFY(this%random)
        CALL rhs_get_random(this%rhs,this%random,stat_info_sub)

        RETURN          
        
      END SUBROUTINE particles_init

      
      SUBROUTINE particles_display_parameters(this,flag,stat_info)
        
        TYPE(Particles), INTENT(IN)     :: this
        INTEGER, INTENT(IN)             :: flag
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: num_species
        INTEGER                         :: rank
        INTEGER                         :: stat_info_sub
        
        
        stat_info     = 0
        
        num_species = physics_get_num_species(this%phys,stat_info_sub)
        rank        = technique_get_rank(this%tech,stat_info_sub)
        
        PRINT *, '------------------Start------------------'
        PRINT *, '     Particles parameters'
        PRINT *, '-----------------------------------------'
        
        IF( flag == 0 ) THEN
           
           PRINT *, "before decomposing..."
           
           WRITE (UNIT=*, FMT=100) &
                "rank, num_part_real", ": ", rank, this%num_part_real
           WRITE (UNIT=*, FMT=100) &
                "rank, num_part_all", ": ", rank, this%num_part_all
           WRITE (UNIT=*, FMT=100) &
                "rank, num_part_ghost", ": ", rank, this%num_part_ghost
           WRITE (UNIT=*, FMT=100) &
                "rank, num_part_fluid", ": ", rank, this%num_part_fluid
           WRITE (UNIT=*, FMT=100) &
                "rank, num_part_wall_solid", ": ", rank, this%num_part_wall_solid
           
        ELSE
           
           PRINT *, "after decomposing..."
           
           WRITE (UNIT=*, FMT=100) &
                "rank, num_part_real", ": ", rank, this%num_part_real
           WRITE (UNIT=*, FMT=100) &
                "rank, num_part_all", ": ", rank, this%num_part_all
           WRITE (UNIT=*, FMT=100) &
                "rank, num_part_ghost", ": ", rank, this%num_part_ghost    
           
        END IF
        
        WRITE (UNIT=*, FMT=100) &
             "rank, num_part_colloid", ": ", rank, this%num_part_colloid
        
        PRINT *, '-------------------End-------------------'
        
100     FORMAT (A25,A3,I4,I7)
        
        RETURN
        
      END SUBROUTINE particles_display_parameters
      
      
