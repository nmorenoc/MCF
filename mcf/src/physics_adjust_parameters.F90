      SUBROUTINE physics_adjust_parameters(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : physics_adjust_parameters
        !----------------------------------------------------
        !
        ! Purpose     : Adjust physics parameters after
        !               they have values, and resolve
        !               the small conflicts.
        !      
        ! Reference   :
        !
        ! Remark      :
        !              1)This is a rouinte of Class Physics,
        !               therefore, the variables of physics
        !               object can be accessed directly, 
        !               however, some pointer variables
        !               are easily allocated if we call 
        !               "_get_" routines.
     
        !
        ! Revisions   : V0.1 23.07.2009, original version.
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
        LOGICAL                         :: read_external
        LOGICAL                         :: Brownian
        INTEGER                         :: cc_lub_type
        INTEGER                         :: cc_repul_type
        INTEGER                         :: cw_repul_type
        INTEGER                         :: cw_lub_type
        INTEGER                         :: adaptive_dt
        INTEGER                         :: num_species
        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(:), POINTER :: min_phys
        REAL(MK), DIMENSION(:), POINTER :: max_phys
        INTEGER                         :: wall_rho_type
        REAL(MK)                        :: wall_layer
        REAL(MK), DIMENSION(:), POINTER :: min_phys_t
        REAL(MK), DIMENSION(:), POINTER :: max_phys_t      
        REAL(MK), DIMENSION(:), POINTER :: dx
        INTEGER                         :: lattice_type
        INTEGER, DIMENSION(:), POINTER  :: num_part_dim
        INTEGER, DIMENSION(:), POINTER  :: num_part_dim_t
        INTEGER                         :: num_part_tot
        INTEGER, DIMENSION(6)           :: num_part_wall
        REAL(MK)                        :: cut_off
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        INTEGER                         :: num_colloid
        
        INTEGER                         :: i
        REAL(MK)                        :: coff
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        NULLIFY(min_phys)
        NULLIFY(max_phys)
        NULLIFY(min_phys_t)
        NULLIFY(max_phys_t)      
        NULLIFY(dx)
        NULLIFY(num_part_dim)
        NULLIFY(num_part_dim_t)
        NULLIFY(bcdef)
        
        num_part_wall(1:6) = 0
        
        coff = 100.0_MK        
        
        !----------------------------------------------------
        ! Get min_phys_t, max_phys_t values from
        ! min_phys and max_phys as basics.        
        !----------------------------------------------------
        
        relax_run    = &
             control_get_relax_run(this%ctrl,stat_info_sub)
        read_external= &
             control_get_read_external(this%ctrl,stat_info_sub)
        Brownian     = &
             control_get_Brownian(this%ctrl,stat_info_sub)
        adaptive_dt  = &
             control_get_adaptive_dt(this%ctrl,stat_info_sub)
       
        num_species  = this%num_species
        num_dim      = this%num_dim
        CALL physics_get_min_phys(this,min_phys,stat_info_sub)
        CALL physics_get_max_phys(this,max_phys,stat_info_sub)
        CALL physics_get_min_phys(this,min_phys_t,stat_info_sub)
        CALL physics_get_max_phys(this,max_phys_t,stat_info_sub) 
        lattice_type = this%lattice_type
        CALL physics_get_num_part_dim(this,num_part_dim, stat_info_sub)
        CALL physics_get_num_part_dim(this,num_part_dim_t, stat_info_sub)
        CALL physics_get_dx(this,dx,stat_info_sub)
        cut_off      = this%cut_off
        CALL physics_get_bcdef(this,bcdef,stat_info_sub)
        
        wall_rho_type = &
             boundary_get_wall_rho_type(this%boundary,stat_info_sub)
        
        !----------------------------------------------------
        ! If one species, should be no colloid.
        !----------------------------------------------------
        
        IF ( num_species == 1 ) THEN
           
           this%num_colloid = 0
           
        END IF
        
        num_colloid  = this%num_colloid
        
        !----------------------------------------------------
        ! Adjust according to different lattice type.
        !----------------------------------------------------
        
        SELECT CASE( lattice_type )
           
           !-------------------------------------------------
           ! If we are using square 2D or simple cubic 3D,
           ! the numbers of particles in each dimension are 
           ! known, then we need to calculate dx.
           ! 
           ! For 2D staggered grid or 3D body center grid
           ! it is the same.
           !-------------------------------------------------
           
        CASE (1:2)
           
           dx(1:num_dim) = &
                (max_phys(1:num_dim) - min_phys(1:num_dim))/ &
                num_part_dim(1:num_dim)
           
           CALL physics_set_dx(this,dx(1:num_dim),stat_info_sub)
           
        CASE (3)
           
           !-------------------------------------------------
           ! 2D Hexagonal lattice.
           !------------------------------------------------
           
           IF ( num_dim == 2 ) THEN
              
              dx(2) = &
                   (max_phys(2) - min_phys(2)) / num_part_dim(2)
              
              dx(1) = dx(2) * SQRT(3.0_MK)
              
              !num_part_dim(1) = &
              !     CEILING((max_phys(1)-min_phys(1)) / dx(1))
              
              !----------------------------------------------
              ! Trying to shrink a little bit in x direction,
              ! in order to make sure particles are evenly 
              ! distributed.
              !----------------------------------------------
              
              !dx(1) = &
              !     (max_phys(1)-min_phys(1))/num_part_dim(1)
              
              CALL physics_set_dx(this,dx(1:num_dim),stat_info_sub)
              
              !----------------------------------------------
              ! Re-calculate the number of particles in each
              ! direction in case we have wall,
              ! take "Ceiling" of the integer, which
              ! might be bigger than the real particles
              ! number.
              !----------------------------------------------
              
              num_part_dim(1:2) = CEILING ( &
                   (max_phys_t(1:2) - min_phys_t(1:2)) / &
                   dx(1:2) )
              
              num_part_dim_t(1:2) = num_part_dim(1:2)

              !----------------------------------------------
              ! 3D face centered lattice.
              !----------------------------------------------
              
           ELSE IF(num_dim == 3) THEN
              
           END IF
           
        CASE DEFAULT
           
           PRINT *, "physics_adjust_parameters : ", &
                "lattice type not available !"
           stat_info = -1
           GOTO 9999
           
        END SELECT ! lattice_type
        
        
        !-------------------------------------------------
        ! 1 Calculate the nubmer of particles in solid
        !   wall, which is handeld by mcf, not by ppm.
        !
        ! 2 Extend the total physical domain to
        !   include the wall particles.
        !
        ! 3 Add up the wall particles to total
        !   particles in each dimension.
        !-------------------------------------------------
        
        DO i = 1, num_dim
           
           IF ( bcdef(2*i-1) == ppm_param_bcdef_wall_solid ) THEN
              
              wall_layer = (1.0_MK+mcf_wall_in_layer_coeff) * cut_off
              
              IF ( wall_rho_type == mcf_wall_rho_type_dynamic ) THEN
                 
                 wall_layer = wall_layer * 2.0_MK
                 
              END IF
              
              num_part_wall(2*i-1) = &
                   CEILING(wall_layer / dx(i))
              
              min_phys_t(i) = min_phys_t(i) - &
                   num_part_wall(2*i-1) * dx(i)
              
              num_part_dim_t(i) = num_part_dim_t(i) + &
                   num_part_wall(2*i-1)
              
           END IF ! bcdef(2i-1)
           
           IF ( bcdef(2*i) == ppm_param_bcdef_wall_solid ) THEN
              
              wall_layer = (1.0_MK+mcf_wall_in_layer_coeff) * cut_off
              
              IF ( wall_rho_type == mcf_wall_rho_type_dynamic ) THEN
                 
                 wall_layer = wall_layer * 2.0_MK
                 
              END IF
              
              num_part_wall(2*i) = &
                   CEILING(wall_layer / dx(i))
              
              max_phys_t(i) = max_phys_t(i) + &
                   num_part_wall(2*i) * dx(i)
              
              num_part_dim_t(i) = num_part_dim_t(i) + &
                   num_part_wall(2*i)
              
           END IF ! bcdef(2i)
           
        END DO ! i = 1, num_dim
           
        !----------------------------------------------------
        ! Reset the min_phys_t, max_phys_t.
        ! Reset the number of particles in each dimension.
        ! Reset the number of particles in total.
        !----------------------------------------------------
        
        CALL physics_set_min_phys_t(this,min_phys_t,stat_info_sub)
        CALL physics_set_max_phys_t(this,max_phys_t,stat_info_sub)
        num_part_tot = PRODUCT(num_part_dim_t(:),1)
        
        CALL physics_set_num_part_dim(this,num_part_dim, stat_info_sub)
        CALL physics_set_num_part_dim_t(this,num_part_dim_t, stat_info_sub)
        CALL physics_set_num_part_tot(this,num_part_tot, stat_info_sub)
        
        CALL boundary_set_min_phys_t(this%boundary,min_phys_t,stat_info_sub)
        CALL boundary_set_max_phys_t(this%boundary,max_phys_t,stat_info_sub)
        
        !----------------------------------------------------
        ! Set the total physics boundary for colloid.
        ! Set the minimal distance of a fluid particle 
        ! from surface of colloids, if there is any.
        !----------------------------------------------------
           
        IF ( num_colloid > 0 ) THEN
           
           cc_lub_type = &
                control_get_cc_lub_type(this%ctrl,stat_info_sub) 
           cc_repul_type = &
                control_get_cc_repul_type(this%ctrl,stat_info_sub)
           cw_lub_type = &
                control_get_cw_lub_type(this%ctrl,stat_info_sub)
           cw_repul_type = &
                control_get_cw_repul_type(this%ctrl,stat_info_sub)
           CALL colloid_set_cc_lub_type(this%colloids,&
                cc_lub_type,stat_info_sub)
           CALL colloid_set_cc_repul_type(this%colloids,&
                cc_repul_type,stat_info_sub)          
           CALL colloid_set_cw_lub_type(this%colloids,&
                cw_lub_type,stat_info_sub)
           CALL colloid_set_cw_repul_type(this%colloids,&
                cw_repul_type,stat_info_sub)
            CALL colloid_set_min_phys_t(this%colloids, &
                min_phys_t,stat_info_sub)
           CALL colloid_set_max_phys_t(this%colloids, &
                max_phys_t,stat_info_sub)
           CALL colloid_set_dout(this%colloids, &
                dx(1)*mcf_colloid_out_layer_coeff,stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Get absolute value of body force and record it.
        !----------------------------------------------------
        
        this%fa_max = SQRT(DOT_PRODUCT(&
             this%body_force(1:num_dim),this%body_force(1:num_dim)))
    
        
        !----------------------------------------------------
        ! Initialize time step.
        !----------------------------------------------------
        
        CALL physics_initialize_dt(this,stat_info_sub)
        
        !-------------------------------------------------
        ! IF step is given from input,
        ! set time end to be negative for adaptive dt.
        !-------------------------------------------------
        
        IF( this%step_end >= 0 ) THEN
           
           IF ( adaptive_dt > 0) THEN
           
              this%time_start = 0.0_MK
              this%time_end   = -1.0_MK
              
           ELSE
              
              this%time_start = this%step_start * this%dt
              this%time_end   = this%step_end * this%dt
              
           END IF
           
           !-------------------------------------------------
           ! IF time is given from input,
           ! set step end to negative for adaptive dt.
           !-------------------------------------------------
           
        ELSE IF ( this%time_end >=0.0_MK ) THEN
              
           IF ( adaptive_dt > 0) THEN
              
              this%step_start = 0
              this%step_end   = -1.0
              
           ELSE
              
              this%step_start = INT(this%time_start / this%dt)
              this%step_end   = CEILING(this%time_end / this%dt)
              
           END IF
         
        END IF ! step_end >= 0.0
           
        
        !----------------------------------------------------
        ! For relax run type 1, set step_relax and time_relax.
        !----------------------------------------------------
        
        IF ( relax_run ) THEN
           
           SELECT CASE ( this%relax_type )
              
           CASE (1)
              
              IF ( this%time_relax > 0 ) THEN
                 
                 this%step_relax = CEILING(this%time_relax / this%dt_relax)
                 
              ELSE IF ( this%step_relax > 0 ) THEN
                 
                 this%time_relax = this%step_relax * this%dt_relax
                 
              END IF
              
              this%disorder_level = 1.0_MK
              
           CASE (2)
              
              this%step_relax = -1
              this%time_relax = -1.0_MK
              
           END SELECT
        
        END IF
        
        !----------------------------------------------------
        ! For non-Brownian simulation, set kt = 0.
        !----------------------------------------------------
        
        IF ( .NOT. Brownian ) THEN
           
           this%kt = 0.0_MK
           
        END IF
        
        !----------------------------------------------------
        ! In case there is periodic or Lees-Edwards
        ! boundary condition, we have to generate the
        ! images of each colloid.
        !----------------------------------------------------
        
        IF ( num_colloid > 0 ) THEN
           
           CALL colloid_adjust_parameters(this%colloids,stat_info_sub)
           CALL colloid_initialize_image(this%colloids,stat_info_sub)
           CALL colloid_compute_image(this%colloids,stat_info_sub)

        END IF
        
9999    CONTINUE       
        
        IF(ASSOCIATED(min_phys))THEN
           DEALLOCATE(min_phys)
        END IF
        
        IF(ASSOCIATED(max_phys))THEN
           DEALLOCATE(max_phys)
        END IF
        
        IF(ASSOCIATED(min_phys_t))THEN
           DEALLOCATE(min_phys_t)
        END IF
        
        IF(ASSOCIATED(max_phys_t))THEN
           DEALLOCATE(max_phys_t)
        END IF
        
        IF(ASSOCIATED(dx))THEN
           DEALLOCATE(dx)
        END IF
        
        IF(ASSOCIATED(num_part_dim)) THEN
           DEALLOCATE(num_part_dim)
        END IF
        
        IF(ASSOCIATED(num_part_dim_t)) THEN
           DEALLOCATE(num_part_dim_t)
        END IF
        
        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF
        
        
        RETURN
        
      END SUBROUTINE physics_adjust_parameters
      
      
