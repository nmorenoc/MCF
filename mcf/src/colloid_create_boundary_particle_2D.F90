      SUBROUTINE colloid_create_boundary_particle_2D(this,&
           dx,p_x,p_sid,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_create_boundary_particle
        !----------------------------------------------------
        !
        ! Purpose     : create the particles on the surface
        !               of colloids;
        !                              
        !
        ! Reference   :
        !
        ! Remark      : 1)
        !               In case of periodic or Lees-Edwards
        !               boundaries, the images of colloid's
        !               center has to be taken into account
        !               to decide if a potential boundary
        !               particle is inside the geometry of
        !               a colloid.
        !               For 2D, maximum 3**2=9=8 images;
        !               For 3D, maximum 3**3=27=26 images
        !               (including itself).
        !
        !               2)
        !               Velocity of boundary particle will
        !               be set to zero, no matter if the
        !               colloid's velocity is not zero,
        !               since it needs to be zero for
        !               relax runs. 
        !               If non-zero velocity(according to
        !               colloid velocity) is needed,
        !               it will be set again after relax
        !               run.
        !
        !               3)
        !               Temporary implementation only disk,
        !               boundary particles mass is always
        !               supposed to be the same as solvent.
        !
        ! Revisions   : V0.3 9.3.2011, implemented model 4,5
        !               i.e., psrm, psfdrm.
        !
        !               V0.2 8.10 2010, set boundary particle
        !               velocity zero, to be able to run
        !               relax runs.
        !
        !               V0.1 11.12 2009, original version.
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
        ! Arguments :
        !
        ! Input
        !
        ! this      : object of colloid class.
        ! dx        : initial distance between particles
        ! p_x       : position.
        ! p_sid     : species ID.
        ! stat_info : status of this routine.
        !----------------------------------------------------

        TYPE(colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: dx
        REAL(MK), DIMENSION(:,:), POINTER       :: p_x
        INTEGER, DIMENSION(:), POINTER          :: p_sid
        INTEGER, INTENT(OUT)                    :: stat_info
        
        
        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------
      
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num_dim
        REAL(MK), DIMENSION(3)                  :: length
        INTEGER                                 :: num_colloid
        INTEGER                                 :: num_layer
        REAL(MK)                                :: rad
        REAL(MK)                                :: angle
        INTEGER                                 :: num_sur
        INTEGER                                 :: num_max
        INTEGER                                 :: num_one
        INTEGER                                 :: num_total

        INTEGER                                 :: i,j,k
        LOGICAL                                 :: counted_boundary
        INTEGER                                 :: shift_layer
        REAL(MK), DIMENSION(3)                  :: sx
        
        REAL(MK), DIMENSION(:,:), POINTER       :: t_x
        INTEGER, DIMENSION(:), POINTER          :: t_sid
        
        REAL(MK), DIMENSION(:,:), POINTER       :: pt_x
        INTEGER, DIMENSION(:), POINTER          :: pt_sid
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        IF ( this%place == mcf_colloid_place_lattice ) THEN
           PRINT *, "colloid_create_boundary_particles_2D: ", &
                "Boundary particles located on lattice !"
           GOTO 9999
        END IF
        
        NULLIFY(t_x)
        NULLIFY(t_sid)

        NULLIFY(pt_x)
        NULLIFY(pt_sid)

        
        !----------------------------------------------------
        ! Get parameters
        !
        ! num_layer : number of layers of boundary particles
        !             around cutoff thickness.
        !----------------------------------------------------
        
        num_dim = this%num_dim
        
        IF ( num_dim /= 2 ) THEN
           PRINT *, "colloid_create_boundary_particle_2D: ", &
                "Dimension should be 2 !"
           stat_info = -1
           GOTO 9999
        END IF
        
        num_colloid = this%num_colloid
        
        num_layer = CEILING(this%din / dx(1))
        
        !PRINT *, "num_layer: ", num_layer
        
        length(1:num_dim) = this%max_phys(1:num_dim) - &
             this%min_phys(1:num_dim)
        
        !----------------------------------------------------
        ! Loop over all colloids and find biggest radius.
        !----------------------------------------------------
        
        rad = this%radius(1,1)
        
        DO j = 2, num_colloid
           
           IF ( this%radius(1,j) > rad ) THEN
              
              rad = this%radius(1,j)
              
           END IF

        END DO
        
        !----------------------------------------------------
        ! Estimate what would be the biggest number of
        ! boundary particles for all colloids.
        ! Each boundary particle's mass is the same and
        ! equal to rho*dx(1)*dx(2)
        !
        ! A:
        ! 1) Suppose we create outest layer with fixed
        ! distance dx(2), x2 for brevity.
        ! 2) Using anti-cosine principle, we get angle
        ! between two adjacent particles at outest layer.
        ! angle1  = ACOS(1.0_MK- (dx(2)/rad)**2/2.0_MK),
        ! 3) From the angle1 we can know the number of particles
        ! num1  on this outest layer.
        !
        ! 4) However, if we suppose the arc length is x2,
        !  then angle2=x2/R will be smaller than angle1,
        ! therefore num2=2*pi/angle2=2*pi*R/x2 > num1
        !
        ! we use num2 as estimated maximal number for each layer,
        ! i.e., num_sur1.
        !
        ! B:
        ! 1) Suppose we create outest layer with total mass
        ! of particles as the shell mass rho*(pi*R^2-pi(R-x1)^2),
        ! i.e., rho*(2*pi*R*x1-pi*x1**2), then the number on this
        ! layer num_sur2 would be (2*pi*R/x2-pi*x1/x2)
        ! 
        ! Since num_sur2 < num_sur1, therefore num_sur=num_sur1
        !
        ! num_max : estimated maximal number of total paritcles..
        !----------------------------------------------------
        
        num_sur = CEILING(2.0_MK*mcf_pi*rad/ dx(2))
        
        num_max = num_colloid * num_layer * num_sur
        
        ALLOCATE(pt_x(num_dim,num_max))
        ALLOCATE(pt_sid(num_max))
        
        !----------------------------------------------------
        ! Set temporary total number of boundary particles 
        ! to zero.
        !----------------------------------------------------
        
        num_total = 0
        
        !----------------------------------------------------
        ! Loop over all colloids
        !----------------------------------------------------
        
        DO j = 1, num_colloid
           
           rad = this%radius(1,j)
           
           !-------------------------------------------------
           ! Estimate what would be the biggest number of
           ! boundary particles for this colloid.
           ! num_max : estimated maximal number.
           ! Then allocate memory for this colloid.
           !-------------------------------------------------
           
           num_sur   = &
                CEILING(2.0_MK*mcf_pi*rad/dx(2))
           
           num_max = num_layer * num_sur
           
           IF (ASSOCIATEd(t_x) ) THEN
              DEALLOCATE(t_x) 
           END IF
           IF (ASSOCIATEd(t_sid) ) THEN
              DEALLOCATE(t_sid) 
           END IF
           ALLOCATE(t_x(num_dim,num_max))
           ALLOCATE(t_sid(num_max))
           
           !-------------------------------------------------
           ! Set number of boundary particles of this 
           ! particular one colloid to zero.
           !-------------------------------------------------
           
           num_one = 0
           
           !-------------------------------------------------
           ! Check its shape
           !-------------------------------------------------
           
           SELECT CASE( this%shape(j) )
              
              !----------------------------------------------
              ! 2D Disk. ( currently only available)
              !----------------------------------------------
              
           CASE (mcf_colloid_shape_cylinder:mcf_colloid_shape_disk)
              
              !----------------------------------------------
              ! Check if placing boundary particles with
              ! fixed distance of each layer , 
              ! fixed number of each layer, or
              ! same mass of each shell.
              !----------------------------------------------
              
              SELECT CASE ( this%place ) 
                 
              CASE (mcf_colloid_place_psfd)
                 
                 !-------------------------------------------
                 ! each layer has the same distance dx(1)
                 ! from adjacent layers.
                 !-------------------------------------------
                 
                 shift_layer  = 0
                 
                 DO i = 0, num_layer-1
                    
                    !----------------------------------------
                    ! Calculate the radius of current layer.
                    !----------------------------------------
                    
                    rad = this%radius(1,j) - i * dx(1)
                    
                    IF ( rad >= dx(2)/2.0_MK ) THEN
                       
                       !-------------------------------------
                       ! Diameter of this layer is >= dx(2)
                       ! Use arccosine theorem to calculate
                       ! angle, whichs give particles 
                       ! distance as dx on surface.
                       !-------------------------------------
                       
                       angle   = &
                            ACOS(1.0_MK- (dx(2)/rad)**2/2.0_MK)
                       num_sur = &
                            CEILING(2.0_MK*mcf_pi/angle)
                       
                       !-------------------------------------
                       ! Adjust angles uniformly.
                       !-------------------------------------
                       
                       angle = 2.0_MK*mcf_pi / num_sur
                       
                       !-------------------------------------
                       ! Loop num_sur position on this layer,
                       ! give a shift of half angle for each
                       ! new layer.
                       !-------------------------------------
                       
                       DO k = 0, num_sur-1
                          
                          sx(1) = this%x(1,j) + rad * &
                               COS(k*angle + 0.5_MK*shift_layer*angle)
                          sx(2) = this%x(2,j) + rad * &
                               SIN(k*angle + 0.5_MK*shift_layer*angle)
                          
                          num_one          = num_one + 1
                          t_x(1:2,num_one) = sx(1:2)
                          t_sid(num_one)   = j
                          
                       END DO ! k = 0, num_sur - 1
                       
                       shift_layer = MOD(shift_layer+1,2)
                       
                       !PRINT *, "num_sur: ", num_sur
                       
                    ELSE
                       
                       !-------------------------------------
                       ! Otherwise, rad is too small already,
                       ! record the center of colloid.
                       !-------------------------------------
                       
                       num_one          = num_one + 1
                       t_x(1:2,num_one) = this%x(1:2,j)
                       t_sid(num_one)   = j
                       
                    END IF ! rad > dx(2)/2.0
                    
                 END DO ! i = 0, num_layer - 1
                 
              CASE (mcf_colloid_place_psfn)
                 
                 !-------------------------------------------
                 ! Each layer has distance dx(1) from adjacent
                 ! layers and
                 ! each layer has the same number of particles.
                 !
                 ! Use arccosine theorem to calculate angle,
                 ! which give particles distance as dx on 
                 ! surface.
                 !-------------------------------------------
                 
                 rad = this%radius(1,j)
                 
                 angle   = &
                      ACOS(1.0_MK-(dx(2)/rad)**2.0_MK/2.0_MK)
                 num_sur = &
                      CEILING(2.0_MK*mcf_pi/angle)
                 
                 !-------------------------------------------
                 ! Adjust angles uniformly.
                 !-------------------------------------------
                 
                 angle = 2.0_MK*mcf_pi / num_sur
                 
                 
                 !-------------------------------------------
                 ! each layer has same distance with dx(2).
                 !-------------------------------------------
                 
                 shift_layer  = 0
                 
                 DO i = 0, num_layer-1
                    
                    rad = this%radius(1,j) - i * dx(1)
                    
                    !----------------------------------------
                    ! Loop num_sur position on this layer, 
                    ! give a shift of half angle for each 
                    ! new layer.
                    !----------------------------------------
                    
                    DO k = 0, num_sur-1
                       
                       sx(1) = this%x(1,j) + rad * &
                            COS(k*angle + 0.5_MK*shift_layer*angle )
                       
                       sx(2) = this%x(2,j) + rad * &
                            SIN(k*angle + 0.5_MK*shift_layer*angle )
                       
                       num_one          = num_one + 1
                       t_x(1:2,num_one) = sx(1:2)
                       t_sid(num_one)   = j
                       
                    END DO ! k = 0, num_sur-1
                    
                    shift_layer = MOD(shift_layer+1,2)
                    
                    !PRINT *, "num_sur: ", num_sur
                    
                 END DO ! i = 0, num_layer
              
                 
              CASE (mcf_colloid_place_psrm)
                 
                 !-------------------------------------------
                 ! each layer has same distance dx(2) from
                 ! adjacent layers.
                 !-------------------------------------------
                 
                 shift_layer  = 0
                 
                 DO i = 0, num_layer-1
                    
                    !----------------------------------------
                    ! Calculate the radius of current layer
                    ! and the mass of this shell.
                    !----------------------------------------
                    
                    rad = this%radius(1,j) - i * dx(1)
                    
                    IF ( rad >= dx(2)/2.0_MK ) THEN
                       
                       !-------------------------------------
                       ! Radius of this layer is >= dx(2)/2.
                       ! 1) rad>=dx(1)
                       ! calculate area of this shell which is
                       ! area = pi*dx(1)*(2.0_MK*rad-dx(1)).
                       ! 2) rad<dx(1)
                       ! calculate area of the disk which is
                       ! area = pi*rad^2
                       ! 3)area unit for a particle is 
                       ! dx(1)*dx(2).
                       ! Then num_sur = area/area_unit
                       ! = pi*(2R-x1)/x2
                       !-------------------------------------
                       
                       IF ( rad >= dx(1) ) THEN
                          num_sur = &
                               CEILING(mcf_pi*(2.0_MK*rad-dx(1))/dx(2))
                       ELSE
                          num_sur = &
                               CEILING(mcf_pi*rad**2/dx(1)/dx(2))
                       END IF
                       
                       !-------------------------------------
                       ! Adjust angles uniformly.
                       !-------------------------------------
                       
                       angle = 2.0_MK*mcf_pi / num_sur
                       
                       !-------------------------------------
                       ! Loop num_sur position on this layer,
                       ! give a shift of half angle for each
                       ! new layer.
                       !-------------------------------------
                       
                       DO k = 0, num_sur-1
                          
                          sx(1) = this%x(1,j) + rad * &
                               COS(k*angle + 0.5_MK*shift_layer*angle)
                          sx(2) = this%x(2,j) + rad * &
                               SIN(k*angle + 0.5_MK*shift_layer*angle)
                          
                          num_one          = num_one + 1
                          t_x(1:2,num_one) = sx(1:2)
                          t_sid(num_one)   = j
                          
                       END DO ! k = 0, num_sur - 1
                       
                       shift_layer = MOD(shift_layer+1,2)
                       
                       !PRINT *, "num_sur: ", num_sur

                    ELSE
                       
                       !-------------------------------------
                       ! Otherwise, rad is too small already,
                       ! record the center of colloid.
                       !-------------------------------------
                       
                       num_one          = num_one + 1
                       t_x(1:2,num_one) = this%x(1:2,j)
                       t_sid(num_one)   = j
                       
                    END IF ! rad > dx(2)/2.0
                    
                 END DO ! i = 0, num_layer - 1
                 
              CASE (mcf_colloid_place_psfdrm)
                 
                 !-------------------------------------------
                 ! each layer has the same distance dx(1)
                 ! from adjacent layers and represent real
                 ! mass of that shell.
                 !-------------------------------------------
                 
                 shift_layer  = 0
                 
                 DO i = 0, num_layer-1
                    
                    !----------------------------------------
                    ! Calculate the radius of current layer.
                    !----------------------------------------
                    
                    rad = this%radius(1,j) - i * dx(1) - dx(1)/2.0_MK
                    
                    IF ( rad >= dx(2)/2.0_MK ) THEN
                       
                       !-------------------------------------
                       ! Diameter of this layer is >= dx(2)
                       ! Use arccosine theorem to calculate
                       ! angle, whichs give particles 
                       ! distance as dx on surface.
                       !-------------------------------------
                       
                       angle   = &
                            ACOS(1.0_MK- (dx(2)/rad)**2/2.0_MK)
                       num_sur = &
                            CEILING(2.0_MK*mcf_pi/angle)
                       
                       !-------------------------------------
                       ! Adjust angles uniformly.
                       !-------------------------------------
                       
                       angle = 2.0_MK*mcf_pi / num_sur
                       
                       !-------------------------------------
                       ! Loop num_sur position on this layer,
                       ! give a shift of half angle for each
                       ! new layer.
                       !-------------------------------------
                       
                       DO k = 0, num_sur-1
                          
                          sx(1) = this%x(1,j) + rad * &
                               COS(k*angle + 0.5_MK*shift_layer*angle)
                          sx(2) = this%x(2,j) + rad * &
                               SIN(k*angle + 0.5_MK*shift_layer*angle)
                          
                          num_one          = num_one + 1
                          t_x(1:2,num_one) = sx(1:2)
                          t_sid(num_one)   = j
                          
                       END DO ! k = 0, num_sur - 1
                       
                       shift_layer = MOD(shift_layer+1,2)
                       
                       !PRINT *, "num_sur: ", num_sur
                       
                    ELSE
                       
                       !-------------------------------------
                       ! Otherwise, rad is too small already,
                       ! record the center of colloid.
                       !-------------------------------------
                       
                       num_one          = num_one + 1
                       t_x(1:2,num_one) = this%x(1:2,j)
                       t_sid(num_one)   = j
                       
                    END IF ! rad > dx(2)/2.0
                    
                 END DO ! i = 0, num_layer - 1
                 
              CASE DEFAULT
                 
                 PRINT *, "colloid_create_boundary_particle_2D : ", &
                      "Type of placing  not supported ! "
                 stat_info = -1
                 GOTO 9999
                 
              END SELECT ! place
              
           CASE DEFAULT
              
              PRINT *, "colloid_create_boundary_particle_2D : ", &
                   "shape not supported ! "
              stat_info = -1
              GOTO 9999
              
           END SELECT ! shape(j)
           
           !-------------------------------------------------
           ! Loop over each generated boundary particle to
           ! check it is inside the domain.
           !-------------------------------------------------
           
           DO i = 1, num_one
              
              counted_boundary = .TRUE.
              
              DO k = 1, num_dim
                 
                 IF ( t_x(k,i) < this%min_phys(k) ) THEN
                    
                    IF ( this%bcdef(2*k-1) == ppm_param_bcdef_periodic .OR. &
                         this%bcdef(2*k-1) == ppm_param_bcdef_LE) THEN
                       
                       t_x(k,i) = t_x(k,i) + length(k)
                       
                    ELSE
                       
                       counted_boundary = .FALSE.
                       
                    END IF ! bcdef
                    
                 END IF ! t_x(k,i) < min_phys.
                    
                 IF ( t_x(k,i) >= this%max_phys(k) ) THEN
                    
                    IF ( this%bcdef(2*k) == ppm_param_bcdef_periodic .OR. &
                         this%bcdef(2*k) == ppm_param_bcdef_LE) THEN
                       
                       t_x(k,i) = t_x(k,i) - length(k)
                       
                    ELSE
                       
                       counted_boundary = .FALSE.
                       
                    END IF ! bcdef
                    
                 END IF !  t_x(k,i) >= max_phys
                 
              END DO ! k
              
              IF ( counted_boundary ) THEN
                 
                 num_total = num_total + 1
                 pt_x(1:num_dim,num_total) = t_x(1:num_dim,i)
                 pt_sid(num_total)         = t_sid(i)
                 
              END IF
              
           END DO ! i = 1, num_one
           
        END DO ! j = 1, num_colloid
        
        !----------------------------------------------------
        ! Allocate memory for output parameters.
        !----------------------------------------------------

        IF ( ASSOCIATED(p_x) ) THEN
           DEALLOCATE(p_x)
        END IF
        IF ( ASSOCIATED(p_sid) ) THEN
           DEALLOCATE(p_sid)
        END IF
        
        ALLOCATE(p_x(num_dim,num_total))
        ALLOCATE(p_sid(num_total))

        !----------------------------------------------------
        ! Copy the result to output parameters.
        !----------------------------------------------------
        
        p_x(1:num_dim,1:num_total) = pt_x(1:num_dim,1:num_total)
        p_sid(1:num_total) = pt_sid(1:num_total)
        
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Release local dynamic memmory.
        !----------------------------------------------------
        
        IF (ASSOCIATEd(t_x) ) THEN
           DEALLOCATE(t_x) 
        END IF

        IF (ASSOCIATEd(t_sid) ) THEN
           DEALLOCATE(t_sid) 
        END IF
        
        IF (ASSOCIATEd(pt_x) ) THEN
           DEALLOCATE(pt_x) 
        END IF

        IF (ASSOCIATEd(pt_sid) ) THEN
           DEALLOCATE(pt_sid) 
        END IF
        
        RETURN
        
      END SUBROUTINE  colloid_create_boundary_particle_2D
      

      
