      SUBROUTINE  particles_init_global_inter_hexagonal(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_init_global_inter_hexagonal
        !----------------------------------------------------
        !
        ! Purpose     : Create positions for particles in 
        !               total physical domain on hexagonal
        !               lattice.
        !              
        !
        ! Routines    :
        !
        ! Remarks     : We generate fluid particle as usual
        !               using hexagonal lattice.
        !               
        ! References  :
        !
        ! Revisions   : V0.1 22.07.2009, original version.
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
    	! Modules :
    	!----------------------------------------------------
        
        USE ppm_module_find_duplicates
        
        !----------------------------------------------------
        ! Arguments :
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER,INTENT(OUT)	                :: stat_info
        
        !----------------------------------------------------
    	! Local variables start here :
        !
        ! num_dim    : number of dimension.
        ! min_phys_t : minimal coordinate of total domain,
        !              including solid wall if there is.
        ! max_phys_t : maximum coordinate of total domain,
        !              including solid wall if there is.
        ! num_part_tot : 
        !              number of estimated physical particles,
        !              if particles are evenly distributed.
        !              Actually numerical particles can be
        !              smaller, since the center of colloids
        !              is usually void without using numerical
        !              particles;
        !              Actually numerical particles can be 
        !              bigger also, since we may generate
        !              boundary particles of colloid on
        !              parallel surfaces.
        ! num_total  : number of actual particles.
        ! num_total_extra :
        !              number of colloid boundary particles
        !              which are paralle layers to surface.
        ! dx         : initial distance between two particles;
        ! cut_off    : cut off of compact support domain.
        !
        ! boundary   : pointer to a boundary obeject.
        ! num_wall_solid : number of solid wall boundaries.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        
        INTEGER                                 :: num_dim
        REAL(MK), DIMENSION(:), POINTER         :: min_phys_t
        REAL(MK), DIMENSION(:), POINTER         :: max_phys_t        
        INTEGER                                 :: num_part_tot
        INTEGER                                 :: num
        REAL(MK), DIMENSION(:), POINTER         :: dx
        REAL(MK), DIMENSION(2)                  :: sx

        
        !----------------------------------------------------
        ! Flags :
        !
        !----------------------------------------------------
        
        INTEGER                                 :: shift_x
        INTEGER                                 :: shift_y
        
	!----------------------------------------------------
    	! Initialization of variables.
    	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        num_dim = this%num_dim
        
        NULLIFY(min_phys_t)
        NULLIFY(max_phys_t)
        NULLIFY(dx)
        
        
        !----------------------------------------------------
        ! Get physics including boundary parameters.
        !----------------------------------------------------
        
        CALL physics_get_min_phys_t(this%phys, &
             min_phys_t,stat_info_sub)
        CALL physics_get_max_phys_t(this%phys, &
             max_phys_t,stat_info_sub)
        num_part_tot = &
             physics_get_num_part_tot(this%phys,stat_info_sub)
        CALL physics_get_dx(this%phys,dx,stat_info_sub)
        
        !----------------------------------------------------
    	! Allocate memory (large enough) for
        ! fluid particles.
        !
	! x   : position
	! v   : velocity
       	! id  : particle ID,  species ID
    	!----------------------------------------------------
        
        num_part_tot = num_part_tot * 2.0_MK

        ALLOCATE(this%x(2,num_part_tot), STAT=stat_info_sub)
        ALLOCATE(this%v(2,num_part_tot), STAT=stat_info_sub)
        ALLOCATE(this%id(this%num_id,num_part_tot), &
             STAT=stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, &
                "particles_init_global_inter_hexagonal : ", &
                "Allocating memory for variables has problem !"
           stat_info = -1
           GOTO 9999
        END IF
	
        
        !----------------------------------------------------
        ! Reset counters to zero.
        !----------------------------------------------------
        
        num = 0
        
        !----------------------------------------------------
        ! Create a new position of particle and check if it 
        ! is a solid wall or colloid boundary particle, 
        ! if there is solid wall or colloid.
        !----------------------------------------------------
     
        sx(2) = min_phys_t(2) + 0.25_MK * dx(2) 
        
        shift_y = 0
        
        DO WHILE( sx(2) < max_phys_t(2) )
           
           sx(1) = min_phys_t(1) + &
                shift_y * 0.5_MK * dx(1) + 0.25_MK * dx(1)
           
           shift_x = 0
           
           DO WHILE( sx(1) < max_phys_t(1) )
              
              num = num + 1
              
              this%x(1:2,num) = sx(1:2)
              
              sx(1) = sx(1) + dx(1)
              
              shift_x = MOD(shift_x + 1, 2)
              
           END DO ! sx(1)
           
           sx(2) = sx(2) + 0.5_MK * dx(2)

           shift_y = MOD(shift_y+1, 2)
           
        END DO ! sx(2)
        
        this%num_part_real  = num
        
        !----------------------------------------------------
        ! Return.
        !----------------------------------------------------
        
9999	CONTINUE	
        
        IF(ASSOCIATED(min_phys_t)) THEN
           DEALLOCATE(min_phys_t) 
        END IF
        
        IF(ASSOCIATED(max_phys_t)) THEN
           DEALLOCATE(max_phys_t) 
        END IF
        
        IF(ASSOCIATED(dx)) THEN
           DEALLOCATE(dx) 
        END IF
        
        RETURN
        
      END SUBROUTINE particles_init_global_inter_hexagonal
      
