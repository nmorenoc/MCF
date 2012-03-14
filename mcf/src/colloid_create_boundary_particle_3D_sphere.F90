      SUBROUTINE colloid_create_boundary_particle_3D_sphere(this,&
           dx,p_x,sid,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_create_boundary_particle_3D_sphere
        !----------------------------------------------------
        !
        ! Purpose     : create the particles on the surface
        !               of a sphere.
        !                              
        !
        ! Reference   :
        !
        ! Remark      : 1)
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
        !
        ! Revisions   : V0.1  Nov. 21, 2011, 
        !               implemented model 5
        !               i.e.,  psfdrm. 
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
        ! sid       : species ID.
        ! stat_info : status of this routine.
        !----------------------------------------------------
        
        TYPE(colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: dx
        REAL(MK), DIMENSION(:,:), POINTER       :: p_x
        INTEGER, INTENT(IN)                     :: sid
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------
      
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num_dim
        INTEGER                                 :: num_layer
        REAL(MK)                                :: rad,rad_max
        REAL(MK)                                :: d_phi, d_theta
        INTEGER                                 :: num_phi, num_theta
        REAL(MK)                                :: phi, theta
        INTEGER                                 :: num_max
        INTEGER                                 :: num
        
        INTEGER                                 :: i,k,m

        
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: t_x
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        IF ( this%place == mcf_colloid_place_lattice ) THEN
           PRINT *, "colloid_create_boundary_particles_3D_sphere: ", &
                "Boundary particles located on lattice !"
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! Get parameters
        !
        ! num_layer : number of layers of boundary particles
        !             around cutoff thickness.
        !----------------------------------------------------
        
        num_dim = this%num_dim
        
        IF ( num_dim /= 3 ) THEN
           PRINT *, "colloid_create_boundary_particle_3D_sphere : ", &
                "Dimension should be 3 !"
           stat_info = -1
           GOTO 9999
        END IF
        
 
        IF ( this%shape(sid) /= mcf_colloid_shape_sphere ) THEN
           PRINT *, "colloid_create_boundary_particle_3D_sphere : ", &
                "shape should be sphere !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Number of layers for boundary particles
        !----------------------------------------------------

        num_layer = CEILING( this%din / dx(3) - 0.5_MK)
        
        !PRINT *, "num_layer: ", num_layer
        
        !----------------------------------------------------
        ! Find biggest radius, i.e., outest layer.
        !----------------------------------------------------
        
        rad_max = this%radius(1,sid)
        
        !----------------------------------------------------
        ! Estimate what would be the biggest number of
        ! boundary particles for one colloid.
        ! Each boundary particle's mass is the same and
        ! equal to rho*dx(1)*dx(2)*dx(3)
        !----------------------------------------------------
        
        d_phi   = dx(1) / rad_max
        num_phi = NINT(mcf_pi / d_phi)+2
        
        d_theta   = dx(2) / rad_max
        num_theta = NINT(mcf_pi / d_theta)+2
        
        num_max =  num_layer * num_phi * num_theta
        
        !----------------------------------------------------
        ! To be safe, we double num_max
        !----------------------------------------------------
       
        num_max = num_max * 2
        
        ALLOCATE(t_x(num_dim,num_max))
        
      
        !----------------------------------------------------
        ! Set temporary total number of boundary particles 
        ! to zero.
        !----------------------------------------------------
        
        num = 0
        
        SELECT CASE ( this%place ) 
           
        CASE ( mcf_colloid_place_psfdrm ) 
                 
           DO m = 1, num_layer
                    
              !----------------------------------------------
              ! The first layer is rad - dx(3)/2 away 
              ! from surface, the second is rad-3*dx(3)/2
              !  ...
              !----------------------------------------------
                    
              rad = this%radius(1,sid) - m*dx(3)+ dx(3)/2.0_MK
                    
              !----------------------------------------------
              ! The first particle is at east pole.
              !----------------------------------------------
                    
              num = num + 1
              
              t_x(1,num) = rad
              t_x(2,num) = 0.0_MK
              t_x(3,num) = 0.0_MK
              
              d_phi   = dx(1) / rad
              num_phi = NINT(mcf_pi / d_phi)
              d_phi   = mcf_pi / DFLOAT(num_phi)
              
              phi = 0.0_MK
              
              DO i = 1, num_phi-1
                 
                 phi       = phi + d_phi
                 d_theta   = dx(2) / ( rad * DSIN(phi))
                 num_theta = NINT(2.0_MK * mcf_pi / d_theta)
                 d_theta   = 2.0_MK * mcf_pi / DFLOAT(num_theta)
                 
                 theta = 0.0_MK
                 
                 DO k = 1, num_theta
                    
                    num = num + 1
                    t_x(1,num) = rad * DCOS(phi)
                    t_x(2,num) = rad * DSIN(phi) * DSIN(theta)
                    t_x(3,num) = rad * DSIN(phi) * DCOS(theta)
                    theta = theta + d_theta
                    
                 END DO ! num_theta
                 
              END DO ! num_phi
              
              !----------------------------------------
              ! The last particle is at west pole
              !----------------------------------------
              
              num = num + 1
              t_x(1,num)= -rad
              t_x(2,num)= 0.0_MK
              t_x(3,num)= 0.0_MK

           END DO ! m = 1, num_layer
           
        CASE DEFAULT
           
           PRINT *, "colloid_create_boundary_particle_3D_sphere: ", &
                "placement not available !"
           stat_info = -1
           GOTO 9999
           
        END SELECT ! place
        
        
        !----------------------------------------------------
        ! Translate each boundary particle with its
        ! colloid center.
        !----------------------------------------------------        
          
        DO i = 1, num_dim
              
           t_x(i,1:num) = t_x(i,1:num) + this%x(i,sid)
           
        END DO
        
        
        !----------------------------------------------------
        ! Allocate memory for output parameters.
        !----------------------------------------------------
        
        IF ( ASSOCIATED(p_x) ) THEN
           DEALLOCATE(p_x)
        END IF
        
        
        ALLOCATE(p_x(num_dim,num))
        
        !----------------------------------------------------
        ! Copy the result to output parameters.
        !----------------------------------------------------
        
        p_x(1:num_dim,1:num) = t_x(1:num_dim,1:num)
        
        
9999    CONTINUE
        
        
        RETURN
        
      END SUBROUTINE  colloid_create_boundary_particle_3D_sphere
      

      
