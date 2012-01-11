      MODULE Class_Physics
        !----------------------------------------------------
      	!  Class      :	Physics
	!----------------------------------------------------
      	!
      	!  Purpose    :
	!> \brief       Variables and corresponding operations
        !>              for physics quantities.
      	!>	   	
        !>              The variable memebers are private.
        !
        !  Remarks    : Since there is a 'SAVE' after 
        !              'IMPLICIT NONE', all procedures using 
        !               this Module share the same copy of 
        !               this Module.
        !
      	!
      	!  References :
     	!
      	!  Revisions  : 0.2 04.03.2010, including relaxation
        !               parameters.
        !
        !               0.1 03.03.2009, original version.
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
   
        
        USE mcf_header
        
        USE Class_Control
        USE Class_Boundary
        USE Class_Colloid
        
          
        IMPLICIT NONE
        SAVE
        

        
        TYPE Physics
           PRIVATE
           
           !-------------------------------------------------
           ! ctrl           : pointer to control obejct.
           !-------------------------------------------------
           
           TYPE(Control), POINTER             :: ctrl
           
           !-------------------------------------------------
           ! num_species    : number of species in simulation.
           !-------------------------------------------------
           
           INTEGER                            :: num_species
           
           !-------------------------------------------------
           ! num_dim        : number of dimension.
           ! min_phys       : minimal boundary of physical domain.
           ! max_phys       : maximal boundary of physical domain.
           ! min_phys_t     : minimal boundary of total physical
           !                  domain, including solid wall, 
           !                  if there is.
           ! max_phys_t     : maximal boundary of total physical
           !                  domain, including solid wall,
           !                  if there is.
           !-------------------------------------------------
           
           INTEGER                            :: num_dim 
           REAL(MK), DIMENSION(:), POINTER    :: min_phys
           REAL(MK), DIMENSION(:), POINTER    :: max_phys
           REAL(MK), DIMENSION(:), POINTER    :: min_phys_t
           REAL(MK), DIMENSION(:), POINTER    :: max_phys_t

     
           !-------------------------------------------------
           ! lattice_type   : lattice type of initial locations
           !                  for particles.
           !                  2D : 
           !                  1=square; 2=staggerd; 3=hexagonal.
           !                  3D : 
           !                  1=simple cubic;
           !                  2=body center;
           !                  3=face center.
           !
           ! num_part_dim   : number of particles per dimension.
           ! num_part_dim_t : number of particles per dimension,
           !                  plus solid wall particles, 
           !                  estimated.
           ! num_part_tot   : number of all particles, 
           !                  estimated.
           !
           ! Note that the aboved three parameters are used
           ! for generating particles internally.
           !
           ! dx             : initial distance between particles.
           ! cut_off        : compact support domain size.
           ! h              : kernel smoothing length.
           !-------------------------------------------------
           
           INTEGER                            :: lattice_type
           INTEGER,  DIMENSION(:), POINTER    :: num_part_dim
           INTEGER,  DIMENSION(:), POINTER    :: num_part_dim_t
           INTEGER                            :: num_part_tot
           REAL(MK), DIMENSION(:), POINTER    :: dx
           REAL(MK)                           :: cut_off
           REAL(MK)                           :: h 
           
           !-------------------------------------------------
           ! dt             : time step size.
           ! dt_c           : time step size by CFL.
           ! dt_nv          : time step size by viscous relaxation.
           ! fa_max         : maximum body force(absolute) per unit mass
           ! dt_f           : time step size by external body force.
           ! step_start     : starting step number.
           ! step_end       : ending step number.
           ! time_start     : starting time.
           ! time_end       : ending time.           
           !-------------------------------------------------
           
           REAL(MK)                           :: dt
           REAL(MK)                           :: dt_c
           REAL(MK)                           :: dt_nu
           REAL(MK)                           :: fa_max
           REAL(MK)                           :: dt_f
           INTEGER                            :: step_start
           INTEGER                            :: step_end
           INTEGER                            :: step_current
           REAL(MK)                           :: time_start
           REAL(MK)                           :: time_end
           REAL(MK)                           :: time_current
           
           !-------------------------------------------------
           ! rho     : initial/required density.
           ! eta     : dynamic/absolute shear viscosity.
           ! ksai    : bulk viscosity.
           ! c       : sound speed.
           ! rho_ref : reference density.
           ! gamma   : power of state equation.
           !-------------------------------------------------

           REAL(MK)                           :: rho
           REAL(MK)                           :: eta
           REAL(MK)                           :: eta_coef
           REAL(MK)                           :: ksai
           REAL(MK)                           :: kt
           REAL(MK)                           :: c
           REAL(MK)                           :: rho_ref
           REAL(MK)                           :: gamma

           !-------------------------------------------------
           ! step_relax :
           ! time_relax :
           ! kt_relax   :
           !-------------------------------------------------
           
           INTEGER                            :: relax_type
           REAL(MK)                           :: dt_relax
           REAL(MK)                           :: dt_c_relax
           INTEGER                            :: step_relax
           REAL(MK)                           :: time_relax
           REAL(MK)                           :: disorder_level
           REAL(MK)                           :: kt_relax
           REAL(MK)                           :: c_relax

           !-------------------------------------------------
           ! Viscoelastic Oldroyd-B model parameters.
           !
           ! tau      : relaxation time of polymer molecules.
	   ! tau_sm   : relaxation time in stochastic model.
           ! n_p      : number of dumbells in one 
           !            SPH/SDPD particle per unit volume.
           ! kt_p     : Boltzmann constant * temperature for
           !            the dumbell.
           ! eigen_dynamics
           !          : indicate if we are using egienvalue &
           !            egienvector dynamics or evolution of
           !            conformation tensor.
           ! eval     : initial egienvalues.
           ! eval_tolerance : if different eigenvalues are different
           !             more than this tolerance, they are treated
           !             as different eigenvalues, otherwise, they
           !             are considered the same.
           !
           ! evec     : initial egienvectors.
           !            2D array notation in order :
           !            ev1_x, ev1_y, ev2_x, ev2_y
           !            3D array notation in order :
           !            ev1_x, ev1_y, ev1_z, ev2_x, ev2_y, ev2_z,
           !            ev3_x, ev3_y, ev3_z
           ! evec_normalize: if eigenvector needed to be normalized
           ! evec_tolerance: if the length of eigenvector exceeds the
           !            unity more than tolerance, it will be normalized.
           !-------------------------------------------------

           REAL(MK)                           :: tau 
           REAL(MK)                           :: n_p
           REAL(MK)                           :: kt_p
           LOGICAL                            :: eigen_dynamics
           REAL(MK), DIMENSION(:), POINTER    :: eval
           REAL(MK)                           :: eval_tolerance
           REAL(MK), DIMENSION(:,:), POINTER  :: evec
           LOGICAL                            :: evec_normalize
           REAL(MK)                           :: evec_tolerance

    	   !-------------------------------------------------
           ! Stochastic modelling
           ! 
	   ! tau_sm   : relaxation time in stochastic model.
           ! k_sm     : turbulence kinetic energy
           !                   
           !-------------------------------------------------

           REAL(MK)                           :: tau_sm
 	   REAL(MK)                           :: k_sm


	   !-------------------------------------------------
           ! body_force_type  : type of body force.
           !                    1 one direction.
           !                    2 two directions shear force.
           !                    3 sin(y) function force.
           !                      e.g., for kolmogorov flow.
	   !			4 stochastic force
           ! body_force       : value of body force.
           ! body_force_d     : differentiaion of body force,
           !                    it is useful when desired flow
           !                    velocity is required.
           !                    It is used to increase or
           !                    decrease every certain number
           !                    of steps(flow_adjust_freq).
           ! flow_direction   : In case we want to have
           !                    constant far field flow 
           !                    velocity, this indicates the
           !                    flow direction.
           ! flow_width       : The width of layer on low 
           !                    boundary of the box and used to
           !                    measure flow velocity.
           ! flow_v           : together with flow_direction
           !                    indicates rquired flow velocity.      
           ! flow_adjust_freq : indicate how frequently to
           !                    change body_force using
           !                    body_force_d to achieve
           !                    required flow velocity.
           !-------------------------------------------------

           INTEGER                            :: body_force_type             
           REAL(MK), DIMENSION(:), POINTER    :: body_force
           REAL(MK), DIMENSION(:), POINTER    :: body_force_d
           INTEGER                            :: flow_direction
           REAL(MK)                           :: flow_width
           REAL(MK)                           :: flow_v
           INTEGER                            :: flow_adjust_freq
           
           !-------------------------------------------------
           ! num_colloid : number of colloid.
           !               If number of species is one,
           !               this will be reset to zero.
           ! colloids    : poiting to a colloid object.
           !-------------------------------------------------
           
           INTEGER                            :: num_colloid
           TYPE(Colloid), POINTER             :: colloids
           
           !-------------------------------------------------
           ! bcdef     : boundary condition definition.
           ! boundary  : pointing a boundary object.
           !-------------------------------------------------
           
           INTEGER, DIMENSION(:), POINTER     :: bcdef
           TYPE(Boundary), POINTER            :: boundary
           
        END TYPE Physics
        
        INTERFACE physics_new
           MODULE PROCEDURE physics_init
        END INTERFACE
        
      CONTAINS
        
#include "physics_new.F90"
#include "physics_finalize.F90"
#include "physics_check_parameters.F90"
#include "physics_adjust_parameters.F90"
#include "physics_initialize_dt.F90"
#include "physics_adapt_dt.F90"
#include "physics_get.F90"
#include "physics_set.F90"
        
        
      END MODULE Class_Physics
      
