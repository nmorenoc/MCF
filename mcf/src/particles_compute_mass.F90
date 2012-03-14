      SUBROUTINE particles_compute_mass(this,stat_info) 
        !----------------------------------------------------
        ! Subroutine  : particles_compute_mass
        !----------------------------------------------------
        !
        ! Purpose     : Computing inital mass of particles,
        !               done on each local process.
        !      
        ! Reference   :
        !
        ! Remark      :
        !               shape = 1: 2D disk / 3D sphere
        !               For mass, it is trival.
        !               For moment of inertia, wikipedia
        !               gives a list.
        !
        !               shape = 2: 2D ellipse / 3D ellipsoid
        !               For moment inertia of ellipse,
        !               see http://www.eformulae.com/
        !               engineering/moment_of_inertia.php
        !               For moment of inertia of ellipsoid,
        !               it is trival, check wikipedia.
        !
        !               shape = 3: 2D star
        !               we calculate its mass and moment
        !               of inertia numerically, i.e.,
        !               take sufficiently thin piece along
        !               radia direction.
        !
        !               shape = 4: 2D / 3D dicolloid, i.e,
        !               two overlapping 2D disks or 3D spheres.
        !               check the report I wrote.
        !
        !
        ! Revisions   : V0.4 21.11.2011, including ellipsoid
        !               and dicollod.
        !
        !               V0.3 17.12.2009, including mass moment
        !               inertia of ellipse.
        !
        !               V0.2 05.10.2009, including mass moment
        !               inertia of disk/sphere.
        !
        !               V0.1 23.07.2009, original version.
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
        ! Arguments :
        !----------------------------------------------------        
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(OUT)		        :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num_dim
        INTEGER                                 :: lattice_type
        REAL(MK), DIMENSION(:), POINTER         :: dx
        REAL(MK), DIMENSION(:), POINTER         :: min_phys
        REAL(MK), DIMENSION(:), POINTER         :: max_phys
        REAL(MK)                                :: init_rho
        REAL(MK)                                :: mass        
        INTEGER, DIMENSION(:), POINTER          :: bcdef
        TYPE(Boundary), POINTER                 :: tboundary
        INTEGER                                 :: num_sym
        INTEGER                                 :: num_colloid
        TYPE(Colloid), POINTER                  :: colloids
        INTEGER, DIMENSION(:), POINTER          :: coll_shape
        REAL(MK), DIMENSION(:,:), POINTER       :: coll_radius
        INTEGER, DIMENSION(:), POINTER          :: coll_freq
        REAL(MK), DIMENSION(:), POINTER         :: coll_m
        REAL(MK), DIMENSION(:,:), POINTER       :: coll_mmi
        REAL(MK)                                :: coll_vol
        REAL(MK)                                :: coll_m_tot
        INTEGER                                 :: i
        REAL(MK)                                :: n_p
        REAL(MK)                                :: d_theta
        REAL(MK)                                :: theta
        REAL(MK)                                :: r
        REAL(MK)                                :: a, b, c, v_cap
        REAL(MK)                                :: ami
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(dx)
        NULLIFY(min_phys)
        NULLIFY(max_phys)
        NULLIFY(bcdef)
        NULLIFY(tboundary)
        NULLIFY(colloids)
        NULLIFY(coll_shape)
        NULLIFY(coll_radius)
        NULLIFY(coll_freq)        
        NULLIFY(coll_m)
        NULLIFY(coll_mmi)
        
        
        !----------------------------------------------------
        ! Dimension and initial density.
        ! Number of symmetry boundaries.
        ! Number of colloid particle.
        !----------------------------------------------------
        
        num_dim = &
             physics_get_num_dim(this%phys,stat_info_sub)
        lattice_type = &
             physics_get_lattice_type(this%phys,stat_info_sub)
        Call physics_get_dx(this%phys,dx,stat_info_sub)
        Call physics_get_min_phys(this%phys,min_phys,stat_info_sub)
        Call physics_get_max_phys(this%phys,max_phys,stat_info_sub)
        init_rho = physics_get_rho(this%phys,stat_info_sub)
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        num_sym = boundary_get_num_sym(tboundary,stat_info_sub)
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        !----------------------------------------------------
        ! Mass of colloids are calculated according
        ! to theire volumens, i.e. mass=v*rho.
        !----------------------------------------------------
        
        coll_m_tot = 0.0_MK
        
        IF ( num_colloid > 0 ) THEN
           
           CAlL physics_get_colloid(this%phys,colloids,stat_info_sub)
           CALL colloid_get_shape(colloids,coll_shape,stat_info_sub)
           CALL colloid_get_radius(colloids,coll_radius,stat_info_sub)
           CALL colloid_get_freq(colloids,coll_freq,stat_info_sub)           
           CALL colloid_get_m(colloids,coll_m,stat_info_sub)
           CALL colloid_get_mmi(colloids,coll_mmi,stat_info_sub)

           coll_vol = 0.0_MK
           
           DO i = 1, num_colloid
              
              IF ( num_dim == 2 ) THEN
                 
                 SELECT CASE(coll_shape(i))
                    
                 CASE( mcf_colloid_shape_cylinder)
                    
                    !----------------------------------------
                    ! 2D cylinder area.
                    !----------------------------------------
                    
                    coll_vol = mcf_pi * coll_radius(1,i)**2.0_MK
                    coll_m(i) = init_rho * coll_vol
                    coll_mmi(3,i) = 0.5_MK*coll_m(i)*coll_radius(1,i)**2
                    
                 CASE( mcf_colloid_shape_disk)
                    
                    !----------------------------------------
                    ! 2D disk area.
                    !----------------------------------------
                    
                    coll_vol = mcf_pi * coll_radius(1,i)**2.0_MK
                    coll_m(i) = init_rho * coll_vol
                    coll_mmi(3,i) = 0.5_MK*coll_m(i)*coll_radius(1,i)**2
             
                    
                 CASE (mcf_colloid_shape_ellipse)
                    
                    !----------------------------------------
                    ! 2D ellipse area.
                    !----------------------------------------
                    
                    coll_vol = mcf_pi*coll_radius(1,i)*coll_radius(2,i)
                    coll_m(i)= init_rho*coll_vol
                    coll_mmi(3,i) = 0.25_MK*coll_m(i)* &
                         ( coll_radius(1,i)**2+coll_radius(2,i)**2 )
                    
                 CASE (mcf_colloid_shape_dicolloid)
                    
                    PRINT *, "particles_compute_mass: ", &
                         "Dicolloid in 2D not implemented."
                    stat_info = -1
                    GOTO 9999
                    
                 CASE (mcf_colloid_shape_star)
                    
                    !----------------------------------------
                    ! 2D star area.
                    !----------------------------------------
                    
                    !coll_vol = mcf_pi * coll_radius(1,i)**2
                    !coll_m(i)= init_rho * coll_vol
                    
                    n_p     = 10000.0_MK
                    d_theta = 2.0_MK*mcf_pi/n_p
                    theta   = d_theta / 2.0_MK
                    
                    coll_vol = 0.0_MK
                    ami      = 0.0_MK
                    
                    DO WHILE ( theta <= 2.0_MK*mcf_pi )
                       
                       r = colloid_polar_star_r(&
                            coll_radius(1,i),coll_radius(2,i), &
                            REAL(coll_freq(i),MK),theta,0.0_MK)
                       
                       coll_vol = coll_vol + r**2 * d_theta / 2.0_MK
                       ami      = ami + r**4 * d_theta / 4.0_MK
                       
                       theta = theta + d_theta
                       
                    END DO
                    
                    coll_m(i)   = coll_vol * init_rho
                    coll_mmi(3,i) = ami * init_rho              
                    
                 CASE  DEFAULT
                    
                    PRINT *, "particles_compute_mass: ", &
                         "No such shape in 2D !"
                    stat_info = -1
                    GOTO 9999
                
                 END SELECT ! shape(i)
                 
                 !-------------------------------------------
                 ! For symmetry boundaries, we reduce the
                 ! mass and moment inertia by factor num_sym.                
                 !-------------------------------------------
                 
                 IF ( num_sym > 0 ) THEN
                    
                    coll_m(i)   = coll_m(i) / REAL(num_sym,MK)
                    coll_mmi(3,i) = coll_mmi(3,i) / REAL(num_sym,MK)
                    
                 END IF
                 
                 
              ELSE IF ( num_dim ==3 ) THEN
                 
                 SELECT CASE(coll_shape(i))
                    
                    !----------------------------------------
                    ! momentum of inertia is not checked!
                    !----------------------------------------
                    
                 CASE(mcf_colloid_shape_cylinder)
                     
                    coll_vol = mcf_pi * &
                         coll_radius(1,i)**2.0_MK * &
                         coll_radius(2,i)
                    coll_m(i) = init_rho * coll_vol
                    coll_mmi(3,i) = 0.5_MK*coll_m(i)*coll_radius(1,i)**2
  
                 CASE(mcf_colloid_shape_sphere)
                    
                    !----------------------------------------
                    ! 3D Sphere volume.
                    !----------------------------------------
                    
                    coll_vol = &
                         4.0_MK* mcf_pi * &
                         coll_radius(1,i)**3.0_MK / 3.0_MK
                    
                    coll_m(i)       = init_rho*coll_vol
                    coll_mmi(1:3,i) = 0.4_MK*coll_m(i)*coll_radius(1,i)**2
                    
                 CASE(mcf_colloid_shape_ellipsoid)
                    
                    !----------------------------------------
                    ! 3D ellipsoid volume.               
                    !----------------------------------------
                    
                    coll_vol = 4.0_MK * mcf_pi * &
                         coll_radius(1,i)*coll_radius(2,i)*&
                         coll_radius(3,i)/3.0_MK
                    
                    coll_m(i) = init_rho * coll_vol
                    coll_mmi(1,i) = 0.2_MK*coll_m(i)* &
                         (coll_radius(2,i)**2+coll_radius(3,i)**2)
                    coll_mmi(2,i) = 0.2_MK*coll_m(i)* &
                         (coll_radius(1,i)**2+coll_radius(3,i)**2)
                    coll_mmi(3,i) = 0.2_MK*coll_m(i)* &
                         (coll_radius(1,i)**2+coll_radius(2,i)**2)
                   
                 CASE(mcf_colloid_shape_dicolloid)
                    !----------------------------------------
                    ! 3D dicolloid, check wikipedia for
                    ! spherical cap to get volume v_cap
                    !----------------------------------------
                    
                    a = coll_radius(1,i)
                    b = coll_radius(2,i)
                    c = SQRT(a**2 - b**2)
                    v_cap = mcf_pi * (2.0_MK * a + b) * &
                         (a - b)**2 / 3.0_MK
                    coll_vol = &
                         8.0_MK * mcf_pi * &
                         a**3.0_MK / 3.0_MK - &
                         2.0_MK * v_cap
                    
                    coll_m(i) = init_rho * coll_vol
                    
                    !----------------------------------------
                    ! check report for details
                    !----------------------------------------
                    
                    coll_mmi(1,i) = &
                         16.0_MK * init_rho * mcf_pi * a**5 / 15.0_MK - &
                         init_rho * mcf_pi * &
                         (a**4*c-2.0_MK*a**2*c**3/3.0_MK+c**5/5.0_MK)
                    
                    coll_mmi(2,i) = init_rho * mcf_pi * &
                         (b**5 + 10.0_MK*a**2*b**3+40.0_MK*a**3*b**2+&
                         45.0_MK*a**4*b+16.0_MK*a**5)/30.0_MK
                    
                    coll_mmi(3,i) = coll_mmi(2,i)
                    
                 CASE  DEFAULT
                    
                    PRINT *, "particles_compute_mass: ", &
                         "No such shape in 3D !"
                    stat_info = -1
                    GOTO 9999
                
                 END SELECT ! coll_shape
                 
                 !-------------------------------------------
                 ! For symmetry boundaries, we reduce the
                 ! mass and moment inertia by factor num_sym.
                 !-------------------------------------------
                 
                 IF ( num_sym > 0 ) THEN
                    
                    coll_m(i)   = coll_m(i) / REAL(num_sym,MK)
                    coll_mmi(1:3,i) = coll_mmi(1:3,i) / REAL(num_sym,MK)
                    
                 END IF
                 
              END IF
              
              coll_m_tot = coll_m_tot + coll_m(i)
              
           END DO ! i = 1 , num_colloid
           
           
           CALL colloid_set_m(colloids,coll_m,stat_info_sub)
           CALL colloid_set_mmi(colloids,coll_mmi,stat_info_sub)
           
        END IF ! num_colloid > 0
        
        IF ( num_dim == 2 ) THEN

           SELECT CASE ( lattice_type )
              
           CASE (mcf_lattice_type_square)
              !----------------------------------------------
              ! Here we consider 2D sqaure or 3D simple
              ! cubic lattice.
              !----------------------------------------------
              
              mass = init_rho
              
              DO i = 1, num_dim
                 
                 mass = mass * dx(i)
                 
              END DO
           
           CASE (mcf_lattice_type_staggered)
              !-------------------------------------------------
              ! 2D staggered lattice.
              !-------------------------------------------------
              
              mass = init_rho
              
              DO i = 1, num_dim
                 
                 mass = mass * dx(i) / SQRT(2.0_MK)
                 
              END DO
              
           CASE (mcf_lattice_type_hexagonal)
              !-------------------------------------------------
              ! 2D hexagonal lattice.
              !-------------------------------------------------
              
              mass = init_rho
              
              DO i = 1, num_dim
                 
                 mass = mass * ( max_phys(i)  - min_phys(i) )
                 
              END DO
              
              mass = mass - coll_m_tot
              
              mass = mass / this%num_part_fluid
              
           CASE DEFAULT
              
              PRINT *, "particles_compute_mass: ",&
                   "lattice type not available !"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! lattice_type
           
        ELSE

           SELECT CASE ( lattice_type )
              
           CASE (mcf_lattice_type_cubic)
              !----------------------------------------------
              ! Here we consider 3D simple cubic lattice.
              !----------------------------------------------
              
              mass = init_rho
              
              DO i = 1, num_dim
                 
                 mass = mass * dx(i)
                 
              END DO
           
           CASE DEFAULT
              
              PRINT *, "particles_compute_mass: ",&
                   "lattice type not available !"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! lattice_type
           
        END IF
        
        this%m(1:this%num_part_real) = mass
        

#if 0   
        
        !------------------------------------------
        ! Alternative way !
        !
        ! Calculate each particle's mass according
        ! to their volumes, i.e. mass=v*rho.
        !------------------------------------------
        
           
        !------------------------------------------
        ! Save the old rhs_density_type, and
        ! set it to 2, i.e. number density
        ! calculation, in order to calculate
        ! particles' mass according to each
        ! volume (1.0/num_density).
        !------------------------------------------
        
        symmetry = &
             control_get_symmetry(this%ctrl,stat_info_sub)
    
        rhs_density_type = &
             control_get_rhs_density_type(this%ctrl,stat_info_sub)
        CALL control_set_rhs_density_type(this%ctrl,2,stat_info_sub)
        CALL rhs_set_rhs_density_type(this%rhs,2,stat_info_sub)
             
        !-----------------------------------------------
        ! For mass density summation 
        ! e.g Morris et al 1997,
        ! it needs other particles' mass, therefore
        ! we have to map ghost of mass.
        !
        ! For number density summation 
        ! e.g. Espanol et al 2003,
        ! it doesn't need other particles' mass, 
        ! therefore we can allocate num_part_real memory.
        !
        ! However, to keep consistency of the interface
        ! of rhs_density_ff(), we map ghost always.
        !-----------------------------------------------        
        
        CALL particles_map_ghost_get(this, &
             l_map_x  = .TRUE., l_map_m = .TRUE., &
             l_map_id = .TRUE.,&
             stat_info=stat_info_sub)
        
        
        !----------------------------------------------------
        ! Get all particles' positions (including ghosts),
      	! to build neighbor list.
      	!----------------------------------------------------
        
        CALL technique_build_list(this%tech,this%x,&
             this%num_part_all,symmetry,stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN
           PRINT *,'particles_compute_mass : ',&
                'Building lists failed !'
           stat_info = -1          
           GOTO 9999           
        END IF
        
        CALL particles_compute_density(this,stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN           
           PRINT *, "particles_compute_mass : ", &
                "Computing number density has problem !"           
           stat_info = -1
           GOTO 9999           
        END IF
        
        IF (symmetry) THEN
           
           CALL particles_map_ghost_put(this, &
                l_map_x = .TRUE., l_map_rho=.TRUE., &
                stat_info=stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "particles_compute_mass : ", &
                   "Receiving number density from ghosts has problem !"
              stat_info = -1
              GOTO 9999           
           END IF
           
        END IF
        
        !--------------------------------
        ! Although the density of ghosts
        ! are not needed, shall we
        ! update ghosts again?
        ! Cause for vgt, it needs, 
        ! otherwise crashed.(12.08.2009)
        !--------------------------------
        
        this%m(1:this%num_part_real) = &
             init_rho / this%rho(1:this%num_part_real)
        
        !------------------------------------------
        ! Restore the old rhs density formulation.
        !------------------------------------------
        
        CALL control_set_rhs_density_type(this%ctrl,&
             rhs_density_type,stat_info_sub)
        CALL rhs_set_rhs_density_type(this%rhs,&
             rhs_density_type,stat_info_sub)
        
#endif
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(dx)) THEN
           DEALLOCATE(dx)
        END IF

        IF(ASSOCIATED(min_phys)) THEN
           DEALLOCATE(min_phys)
        END IF

        IF(ASSOCIATED(max_phys)) THEN
           DEALLOCATE(max_phys)
        END IF
        
        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF
        
        IF(ASSOCIATED(coll_shape)) THEN
           DEALLOCATE(coll_shape)
        END IF
        
        IF(ASSOCIATED(coll_radius)) THEN
           DEALLOCATE(coll_radius)
        END IF
        
        IF(ASSOCIATED(coll_freq)) THEN
           DEALLOCATE(coll_freq)
        END IF
        
        IF(ASSOCIATED(coll_m)) THEN
           DEALLOCATE(coll_m)
        END IF

        IF(ASSOCIATED(coll_mmi)) THEN
           DEALLOCATE(coll_mmi)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_compute_mass
      
      
