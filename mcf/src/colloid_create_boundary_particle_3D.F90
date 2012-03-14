      SUBROUTINE colloid_create_boundary_particle_3D(this,&
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
        !               Temporary implementation only sphere,
        !               boundary particles mass is always
        !               supposed to be the same as solvent.
        !
        ! Revisions   : V0.1 Nov. 21, 2011, implemented model 5
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
        REAL(MK)                                :: rad_max
        REAL(MK)                                :: d_phi, d_theta
        INTEGER                                 :: num_phi, num_theta
        INTEGER                                 :: num_max
        INTEGER                                 :: num_one
        INTEGER                                 :: num_total

        INTEGER                                 :: i,j,k
        LOGICAL                                 :: counted_boundary
        
        REAL(MK), DIMENSION(:,:), POINTER       :: t_x
        
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: pt_x
        INTEGER, DIMENSION(:), ALLOCATABLE      :: pt_sid
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        IF ( this%place == mcf_colloid_place_lattice ) THEN
           PRINT *, "colloid_create_boundary_particles: ", &
                "Boundary particles located on lattice !"
           GOTO 9999
        END IF
        
        NULLIFY(t_x)
     
        !----------------------------------------------------
        ! Get parameters
        !
        ! num_layer : number of layers of boundary particles
        !             around cutoff thickness.
        !----------------------------------------------------
        
        num_dim = this%num_dim
        
        IF ( num_dim /= 3 ) THEN
           PRINT *, "colloid_create_boundary_particle : ", &
                "Dimension should be 3 !"
           stat_info = -1
           GOTO 9999
        END IF
        
        length(1:num_dim) = this%max_phys(1:num_dim) - &
             this%min_phys(1:num_dim)
        
        num_colloid = this%num_colloid
        
        !----------------------------------------------------
        ! Number of layers for boundary particles.
        !----------------------------------------------------

        num_layer = CEILING( this%din / dx(3))
        
        !----------------------------------------------------
        ! Find biggest radius.
        !----------------------------------------------------
        
        rad_max = MAXVAL(this%radius(1,:))
        
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
        
        !----------------------------------------------------
        ! Estimate what would be the biggest number of
        ! boundary particles for all colloids.
        !----------------------------------------------------
    
        num_max = num_colloid * num_layer * num_phi * num_theta
        
        !----------------------------------------------------
        ! To be safe, we triple num_max
        !----------------------------------------------------
        
        num_max = num_max * 3
        
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
           
           SELECT CASE(this%shape(j))
           
           CASE ( mcf_colloid_shape_sphere ) 
              
              CALL colloid_create_boundary_particle_3D_sphere(this,&
                   dx(1:num_dim),t_x,j,stat_info_sub)

              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "colloid_create_boundary_particle_3D: ", &
                      "boundary particles generated for sphere failed !"
                 stat_info = -1
                 GOTO 9999                 
              END IF
              
           CASE ( mcf_colloid_shape_dicolloid ) 

              CALL colloid_create_boundary_particle_3D_dicolloid(this,&
                   dx(1:num_dim),t_x,j,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "colloid_create_boundary_particle_3D: ", &
                      "boundary particles generated for dicolloid failed !"
                 stat_info = -1
                 GOTO 9999                 
              END IF
              
           CASE DEFAULT
              PRINT *, "colloid_create_boundary_particle_3D: ", &
                   "shape not available !"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! shape(j)
           
           num_one = SIZE(t_x,2)
           
           
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
                    
                 ELSE IF ( t_x(k,i) >= this%max_phys(k) ) THEN
                    
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
                 pt_sid(num_total)         = j
                 
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
        p_sid(1:num_total)         = pt_sid(1:num_total)
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Release local dynamic memmory.
        !----------------------------------------------------
        
        IF (ASSOCIATEd(t_x) ) THEN
           DEALLOCATE(t_x) 
        END IF

        RETURN
        
      END SUBROUTINE  colloid_create_boundary_particle_3D
      
#include "colloid_create_boundary_particle_3D_sphere.F90"
#include "colloid_create_boundary_particle_3D_dicolloid.F90"
