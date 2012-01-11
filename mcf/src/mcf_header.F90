      MODULE mcf_header
	!----------------------------------------------------
      	!       mcf_header
      	!----------------------------------------------------
        !
      	!  Purpos     : Parameters for global constants.
      	!
      	!  Remarks    :
      	!
      	!  References :
      	!
      	!  Version    : V 0.1. 03.03.2009, original version.
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
      	! Incude files
        !----------------------------------------------------
        
        IMPLICIT NONE
        SAVE
        
        INCLUDE 'ppm_param.inc'

        
        !----------------------------------------------------
        ! Float number Precision.
        !----------------------------------------------------
        
        !INTEGER, PARAMETER      :: MK        = ppm_kind_single
        INTEGER, PARAMETER      ::  MK        =  ppm_kind_double
        
        !----------------------------------------------------
        ! Default maximum length of string.
        !----------------------------------------------------
        
        INTEGER, PARAMETER      ::  MAX_CHAR   = 256
        
        !------------------------------------------
        ! Pi defined here explicitly.
	!------------------------------------------
        
        REAL(MK), PARAMETER     :: mcf_pi = 3.141592653589793238_MK
        
        !----------------------------------------------------
       	! Zero down to Machine precision
        ! used in mcf is definded here.
        ! epsilon : 2.220446049250313E-016 for double
        ! epsilon : 1.1920929E-07 for single
       	!----------------------------------------------------
        
       	!REAL (MK), PARAMETER	:: mcf_machine_zero = EPSILON(mcf_pi) * 2.0_MK
        REAL (MK)      	        :: mcf_machine_zero

        !----------------------------------------------------
        ! Boltzmann constant defined here,
        ! too small not used yet.
	!----------------------------------------------------
        
        REAL(MK), PARAMETER     :: mcf_kB =1.380650424e-23
        
        !----------------------------------------------------
        ! Maximum absolute value of a random variable.
	!----------------------------------------------------
        
        REAL(MK), PARAMETER     :: mcf_random_Gaussian_max_abs = 6.0_MK

        !----------------------------------------------------
        ! For relax run:
        ! After relax run with Brownian motion, particles
        ! have some potential energy(e.g., some paris are
        ! very close by, which generates big repulsive force
        ! once the real simulation starts.
        ! This potential energy may be 'represented' by the
        ! total kinetic energy of the system. Therefore,
        ! after Brownian motion relaxation, switch T=0 and
        ! run the code further with presure and viscous force
        ! to reduce the kinetic energy until certain 
        ! amount of before.
        !
        ! mcf_kinetics_reduce_factor : indicated factor needs
        ! to be reduced.
        !----------------------------------------------------
        
        REAL(MK), PARAMETER     :: mcf_kinetics_reduce_factor = 100.0_MK
        
        !----------------------------------------------------
        ! Define lattice type
        ! 2D
        ! 1: square 2: staggered 3: hexagonal
        ! 3D
        ! 1: cubic  2: body center 3: face center
        !----------------------------------------------------

        INTEGER, PARAMETER      :: mcf_lattice_type_square    = 1
        INTEGER, PARAMETER      :: mcf_lattice_type_staggered = 2
        INTEGER, PARAMETER      :: mcf_lattice_type_hexagonal = 3

        INTEGER, PARAMETER      :: mcf_lattice_type_cubic     = 1
        INTEGER, PARAMETER      :: mcf_lattice_type_body      = 2
        INTEGER, PARAMETER      :: mcf_lattice_type_face      = 3
        !----------------------------------------------------
        ! Define type of particles
        ! 0         : fluid particle
        ! -2*D - -1 : wall boundary particles
        !             D being number of dimensions.
        ! 1-N       : colloidal boundary particle
        !             N being number of colloids.
        !----------------------------------------------------
        
        INTEGER, PARAMETER      :: mcf_particle_type_fluid   = 0
        INTEGER, PARAMETER      :: mcf_particle_type_wall    = -1
        INTEGER, PARAMETER      :: mcf_particle_type_colloid = 1

        !----------------------------------------------------
        ! Thickness coefficient of inner ring
        ! layer of a colloid,
        ! i.e.
        !    num_part_in_ring = 
        !    ceiling(mcf_colloid_in_layer*cut_off+ cut_off )
        !    inner_ring = num_part_inner_ring * dx
	!----------------------------------------------------
        
        REAL(MK), PARAMETER     :: mcf_colloid_in_layer_coeff = 0.05_MK
      
        !----------------------------------------------------
        ! Thickness coefficient of outside layer of a colloid,
        ! i.e., if a colloid boundary particle exceeds a 
        ! distance of radius bigger than
        ! mcf_colloid_out_layer_coeff * cut_off,
        ! then it is moving inconsistently.
        !----------------------------------------------------
        
        REAL(MK), PARAMETER     :: mcf_colloid_out_layer_coeff = 0.05_MK
  
        !----------------------------------------------------
        ! Maximal ratio of distance between colloid boundary 
        ! particle and fluid particle, for extrapolation, 
        ! i.e.  mcf_colloid_dist_ratio = 
        ! colloid_boundary_particle_surface_dist / 
        ! fluid_particle_surface_dist
        ! 0.5 is the valued used by Morris et al. 1997.
        !----------------------------------------------------
        
        REAL(MK), PARAMETER     :: mcf_colloid_dist_ratio = 0.5_MK
        
      
        
        !----------------------------------------------------
        ! Define the density type of colloid.
        ! 0 constant, no calculation.
        ! 1 two cut_off layers and caculation as solvent
        !----------------------------------------------------

        INTEGER, PARAMETER       :: mcf_colloid_rho_type_constant = 0
        INTEGER, PARAMETER       :: mcf_colloid_rho_type_dynamic  = 1
        
        !----------------------------------------------------
        ! Define the shape of a colloidal particle as
        ! constant, which facilitates programming.
        !
        !----------------------------------------------------

        INTEGER, PARAMETER      :: mcf_colloid_shape_disk      = 1
        INTEGER, PARAMETER      :: mcf_colloid_shape_sphere    = 1
        INTEGER, PARAMETER      :: mcf_colloid_shape_ellipse   = 2
        INTEGER, PARAMETER      :: mcf_colloid_shape_ellipsoid = 2
        INTEGER, PARAMETER      :: mcf_colloid_shape_star      = 3
        INTEGER, PARAMETER      :: mcf_colloid_shape_dicolloid = 4

        !----------------------------------------------------
        ! Define how to place particles for 
        ! a colloidal particle. In four following cases,
        ! each boundary particle has always the same mass
        ! as a solvent particle.
        ! 1: on lattice;
        ! 2: parallel to surface with fixed distance(psfd),
        !    starting from R, then R-dx(1)...
        ! 3: parallel to surface with fixed number of particles
        !    (psfn), starting from R, then R-dx(1)...
        ! 4: parallel to surface with total mass equal to
        !    the real mass on that layer. Therefore, number of
        !    boundary particles is adaptive on each layer,
        !    starting from R, then R-dx(1)...
        ! 5: parallel to surface with total mass equal to
        !    the real mass on that layer. Therefore, number of
        !    boundary particles is adaptive on each layer...
        !    starting from R-dx(1)/2, then R-3dx(1)/2...
        !    it can be considered a combination of 2+4.
        !----------------------------------------------------
        
        INTEGER, PARAMETER      :: mcf_colloid_place_lattice = 1
        INTEGER, PARAMETER      :: mcf_colloid_place_psfd    = 2
        INTEGER, PARAMETER      :: mcf_colloid_place_psfn    = 3
        INTEGER, PARAMETER      :: mcf_colloid_place_psrm    = 4
        INTEGER, PARAMETER      :: mcf_colloid_place_psfdrm  = 5

        !----------------------------------------------------
        ! Define the how two colloidal particles interact.
        ! 0: no interaction at all.
        ! 1: lubrication theory frist order.
        !    ( a colloid - b colloid)
        !----------------------------------------------------
        
        INTEGER, PARAMETER      :: mcf_cc_lub_type_no    = 0
        INTEGER, PARAMETER      :: mcf_cc_lub_type_first = 1
        
         !----------------------------------------------------
        ! Define repulsive force between how two colloidal
        ! particles.
        ! 0: no repulsive at all.
        ! 1: Hookean spring type repulsive force.
        !    ( a colloid - b colloid)
        ! 2: DLVO type repulsive force.
        !----------------------------------------------------
        
        INTEGER, PARAMETER      :: mcf_cc_repul_type_no      = 0
        INTEGER, PARAMETER      :: mcf_cc_repul_type_Hookean = 1
        INTEGER, PARAMETER      :: mcf_cc_repul_type_DLVO    = 2

        !----------------------------------------------------
        ! Define the how colloid interact with wall
        ! 0: no interaction at all.
        ! 1: lubrication theory first order.
        !    ( a colloid - b wall)
        !----------------------------------------------------
        
        INTEGER, PARAMETER      :: mcf_cw_lub_type_no    = 0
        INTEGER, PARAMETER      :: mcf_cw_lub_type_first = 1
        
        !----------------------------------------------------
        ! Define repulsive force between colloidal and wall
        ! particles.
        ! 0: no repulsive at all.
        ! 1: Hookean spring force.
        !    ( a colloid - bwall)
        ! 2: DLVO type repulsive force.
        !----------------------------------------------------
        
        INTEGER, PARAMETER      :: mcf_cw_repul_type_no      = 0
        INTEGER, PARAMETER      :: mcf_cw_repul_type_Hookean = 1
        INTEGER, PARAMETER      :: mcf_cw_repul_type_DLVO    = 2


        !----------------------------------------------------
        ! Thickness coefficient of layer of a wall boundary
        ! i.e. 
        !     num_part_wall(number of particle inside wall)
        !     =  Ceiling( (1+mcf_wall_layer_coeff) * cut_off/ dx )
        !     wall_layer    = num_part_wall * dx
	!----------------------------------------------------
        
        REAL(MK), PARAMETER     :: mcf_wall_in_layer_coeff = 0.05_MK

        !----------------------------------------------------
        ! Maximal ratio of distance between wall boundary 
        ! particle and fluid particle, for extrapolation, 
        ! i.e.  mcf_wall_dist_ratio = 
        !  wall_boundary_particle_surface_dist / 
        ! fluid_particle_surface_dist
        !----------------------------------------------------
        
        !REAL(MK), PARAMETER     :: mcf_wall_dist_ratio = 10.0_MK
        REAL(MK), PARAMETER     :: mcf_wall_dist_ratio = 0.5_MK
        
        !----------------------------------------------------
        ! Define the density type of wall.
        ! 0 constant, no calculation.
        ! 1 two cut_off layers and caculation as solvent
        !----------------------------------------------------

        INTEGER, PARAMETER       :: mcf_wall_rho_type_constant = 0
        INTEGER, PARAMETER       :: mcf_wall_rho_type_dynamic  = 1
        
        !----------------------------------------------------
        ! Indices of two cells interacting with each other,
        ! sym_inter3d : symmetric interaction in 3D.
        ! "0" is the center or first cell.
        ! "13" is the last cell.
        !----------------------------------------------------
        
        INTEGER, DIMENSION(2,14)        :: sym_inter3d = &
             RESHAPE((/0,0,0,1,0,3,0,4,0,9,0,10,0,12,0,13,&
             1,3,1,9,1,12,3,9,3,10,4,9 /),(/2,14/))
        
        !----------------------------------------------------
        ! Relative shift of cell i(from 0 to 13) to 
        ! center/first cell 0.
   	!----------------------------------------------------

        INTEGER, DIMENSION(3,0:13)      :: sym_inp3d  = &
             RESHAPE((/0,0,0,1,0,0,-1,1,0,0,1,0,1,1,0, &
             -1,-1,1, 0,-1,1, 1,-1,1,&
             -1, 0,1, 0, 0,1, 1, 0,1,&
             -1, 1,1, 0, 1,1, 1, 1,1/),(/3,14/))

      END MODULE MCF_HEADER
      
      
