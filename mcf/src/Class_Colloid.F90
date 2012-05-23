      MODULE Class_Colloid
        !----------------------------------------------------
        ! Class	      :	Colloid
        !----------------------------------------------------
        !
        !  Purpose    :
        !> \brief       Variables and corresponding operations
        !>              for colloids quantities.
        !>	   	
        !>              The variable memebers are private.
        !
        !  References :
        !
        !  Remarks    :  In context of particle methods,
        !                there is always a ambuity between
        !                a physical colloid particle and
        !                a colloid boundary particle.
        !                The former is a physical soloid
        !                mesoscopic particle;
        !                The latter is a numerical SPH/DPD/SDPD
        !                particle, denoted as boundary particle,
        !                which is a interpolation point in 
        !                context of SPH, 
        !                is a thermodynamics subsystem in 
        !                context of SDPD.
        !                To distinguish, a colloid always means
        !                the former, i.e. physical particle;
        !                a colloid particle/boundary particle
        !                always mean the latter, i.e. SPH/DPD/SDPD
        !                particle.
        !
        !
        !  Revisions   : V0.2 10.09.2009, add body force(per unit mass)
        !                in Colloid Class, in order to be able
        !                to have seperate body force from fluid.
        !
        !                V0.1 03.03.2009, original version
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
        USE Class_Technique
        USE Class_Boundary
        USE Class_Tool
        
        IMPLICIT NONE
        SAVE
        
        TYPE Colloid

           PRIVATE
           
           TYPE(Technique), POINTER             :: tech

           !-------------------------------------------------
           !  Basic quantities, need to be saved
           !  if restart is required.
           !
           !  Command characteristics:
           !
           !  rho        : density
           !  rho_type   : density type 
           !               0 constant, 1 calculated
           !  translate  : free to translate or not
           !  rotate     : free to rotate or not
           !  place      : how to place boundary particles
           !  noslip_type      : how to model no slip condition
           !  body_force_type  : type of body force.
           !  body_force       : body force per unit mass
           !  cc_lub_cut_off   : cut off surface distance of 
           !                     direct lubrication interaction 
           !                     between colloid and colloid.
           !  cc_lub_cut_on    : cut on surface distance of 
           !                     direct lubrication interaction 
           !                     between colloid and colloid.
           !                     i.e., minimal gap.
           !  cc_repul_cut_off : cut off for repulsive force
           !                     which prevent overlap.
           !  cc_repul_F_max   : maximum repulsive force.
           !  cw_lub_cut_off   : cut off surface distance of 
           !                     direct lubrication interaction 
           !                     between colloid and wall.
           !
           !  h                : smoothing length for caculating
           !                     time step.
           !  dt_f             : time step constranit due to
           !                     maximum acceleration.
           !-------------------------------------------------

           INTEGER                            :: num_dim
           INTEGER                            :: num_colloid
           INTEGER                            :: integrate_type
           INTEGER                            :: integrate_RK
           INTEGER                            :: integrate_AB
           REAL(MK)                           :: adapt_t_coef
           INTEGER                            :: sub_time_step
           INTEGER                            :: implicit_pair_num_sweep
           INTEGER                            :: explicit_sub_time_step
           REAL(MK)                           :: rho
           INTEGER                            :: rho_type
           LOGICAL                            :: translate
           LOGICAL                            :: rotate
           INTEGER                            :: place
           INTEGER                            :: noslip_type
           INTEGER                            :: body_force_type
           REAL(MK), DIMENSION(:), POINTER    :: body_force
           
           INTEGER                            :: cc_lub_type
           REAL(MK)                           :: cc_lub_cut_off
           REAL(MK)                           :: cc_lub_cut_on

           INTEGER                            :: cc_repul_type
           REAL(MK)                           :: cc_repul_cut_off
           REAL(MK)                           :: cc_repul_cut_on
           REAL(MK)                           :: cc_repul_F0
           
           INTEGER                            :: cw_lub_type
           REAL(MK)                           :: cw_lub_cut_off
           REAL(MK)                           :: cw_lub_cut_on

           INTEGER                            :: cw_repul_type
           REAL(MK)                           :: cw_repul_cut_off
           REAL(MK)                           :: cw_repul_cut_on
           REAL(MK)                           :: cw_repul_F0
           
           REAL(MK)                           :: h
           REAL(MK)                           :: dt_f

           !-------------------------------------------------
           !  Individual characteristics:
           !
           !  shape      : shape of colloids
           !  radius     : radius.
           !  freq       : roughness
           !  m          : mass
           !  mmi        : mass momentum of inertia
           !  x          : position
           !  v          : translating velocity
           ! v(:,:,1)    : current time step
           ! v(:,:,2)    : previous time step
           ! v(:,:,3)    : pre-previous time step
           ! v(:,:,4)    : pre-pre-previous time step
           !  drag       : drag
           !  rot_vector : current rotation vector
           !  acc_vector : accumulation rotation vector       
           !  rot_matrix : current rotation matrix
           !  acc_matrix : accumulation rotation matrix
           !  theta      : accumulated rotational angle.
           !               usefull in 2D but not in 3D.
           !  omega      : rotating velocity, 
           !               pseudovector always 3D
           !  about omega(:,:,i), check for v(:,:,i)
           !  torque     : torque, pseudovector always 3D
           !  num_numerical_part : 
           !               number of particles one colloid consist of
           !  m          : mass
           !  mm         : mass momentum inertia
           !-------------------------------------------------

           INTEGER, DIMENSION(:), POINTER     :: shape
           REAL(MK), DIMENSION(:,:), POINTER  :: radius
           INTEGER, DIMENSION(:), POINTER     :: freq
           REAL(MK), DIMENSION(:), POINTER    :: m
           REAL(MK), DIMENSION(:,:), POINTER  :: mmi
           REAL(MK), DIMENSION(:,:), POINTER  :: x
           REAL(MK), DIMENSION(:,:,:),POINTER :: v
           REAL(MK), DIMENSION(:,:), POINTER  :: drag_lub
           REAL(MK), DIMENSION(:,:), POINTER  :: drag_repul
           REAL(MK), DIMENSION(:,:), POINTER  :: drag
           REAL(MK), DIMENSION(:,:), POINTER  :: rot_vector
           REAL(MK), DIMENSION(:,:), POINTER  :: acc_vector           
           REAL(MK), DIMENSION(:,:,:), POINTER:: rot_matrix
           REAL(MK), DIMENSION(:,:,:), POINTER:: acc_matrix
           REAL(MK), DIMENSION(:,:), POINTER  :: theta
           REAL(MK), DIMENSION(:,:,:),POINTER :: omega
           REAL(MK), DIMENSION(:,:), POINTER  :: torque
           INTEGER, DIMENSION(:), POINTER     :: num_physical_part
           INTEGER, DIMENSION(:), POINTER     :: num_numerical_part
           
           !-------------------------------------------------
           !  Derived quantities, can always be
           !  calcualted from basic quantities.
           !  f        : force per unit mass(for translation)
           !  about f(:,:,i), check for v(:,:,i)
           !  fa_min   : mimimum force acceleration.
           !  fa_max   : maximum force acceleration.
           !  alpha    : force per unit mass(for rotation),
           !              pseudovector always 3D
           !  about alpha(:,:,i), check for v(:,:,i)           
           !  k_energy : kinetic energy 
           !  mom      : momentum of  colloids
           !  k_energy_tot : 
           !  mom_tot  :
           !-------------------------------------------------
           
           REAL(MK), DIMENSION(:,:,:),POINTER :: f
           REAL(MK)                           :: fa_min
           REAL(MK)                           :: fa_max
           REAL(MK), DIMENSION(:,:,:),POINTER :: alpha
           REAL(MK), DIMENSION(:), POINTER    :: k_energy
           REAL(MK), DIMENSION(:,:), POINTER  :: mom
           REAL(MK)                           :: k_energy_tot
           REAL(MK), DIMENSION(:), POINTER    :: mom_tot
           INTEGER                            :: num_physical_part_tot
           INTEGER                            :: num_numerical_part_tot
           
           !-------------------------------------------------
           ! Physics parameters :
           !
           ! min_phys : minimal boundary of physical domain.
           ! max_phys : maximal boundary of physical domain.
           ! bcdef    : boundary conditions at physical domain.
           ! boundary : boundary object of Boundary Class.
           ! dout     : minimal distance of a fluid particle from surface.
           ! din      : maximum distance of a boundary particle from surface.
           ! num_image: number of image boxes, according
           !            to different boundary conditions.
           ! x_image : images' positions, including itself.
           ! v_image : images' velocities, including itself.
           !-------------------------------------------------
           
           REAL(MK), DIMENSION(:), POINTER    :: min_phys
           REAL(MK), DIMENSION(:), POINTER    :: max_phys
           REAL(MK), DIMENSION(:), POINTER    :: min_phys_t
           REAL(MK), DIMENSION(:), POINTER    :: max_phys_t
           INTEGER, DIMENSION(:), POINTER     :: bcdef
           TYPE(Boundary), POINTER            :: boundary

           REAL(MK)                           :: cut_off
           REAL(MK)                           :: dout
           REAL(MK)                           :: din
           REAL(MK)                           :: eta

           INTEGER                            :: num_image
           REAL(MK), DIMENSION(:,:,:), POINTER:: x_image
           REAL(MK), DIMENSION(:,:,:), POINTER:: v_image
           
           TYPE(Tool)                         :: tool
           
           
        END TYPE Colloid
        
        INTERFACE colloid_new
           MODULE PROCEDURE colloid_init_default
           MODULE PROCEDURE colloid_init
        END INTERFACE
                
      CONTAINS       
        
        REAL(MK) FUNCTION colloid_polar_angle(x,y)
          !--------------------------------------------------
          ! Get angle[0:2pi) between point (x,y) and 
          ! x+ direction.
          ! ACOS returns a value in [0:pi], therefore,
          ! if y coordinate is negative, we do 2*pi-angle.
          !--------------------------------------------------
          
          REAL(MK),INTENT(IN)           :: x
          REAL(MK),INTENT(IN)           :: y
          
          !--------------------------------------------------
          ! If x=0 && y=0, return angle=0
          !--------------------------------------------------
          
          IF ( ABS(x) <=mcf_machine_zero .AND. &
               ABS(y) <=mcf_machine_zero ) THEN
             
             colloid_polar_angle = 0.0_MK
             
          ELSE
             
             colloid_polar_angle = ACOS(x/SQRT(x**2+y**2))
             
             IF ( y < 0.0_MK ) THEN
                colloid_polar_angle = &
                     2.0_MK*mcf_pi - colloid_polar_angle
             END IF
             
          END IF
          
        END FUNCTION colloid_polar_angle
        
#include "colloid_new.F90"
#include "colloid_finalize.F90"
#include "colloid_check_parameters.F90"
#include "colloid_adjust_parameters.F90"
#include "colloid_get.F90"
#include "colloid_set.F90"
#include "colloid_check_boundary_particle.F90"
#include "colloid_create_boundary_particle_2D.F90"
#include "colloid_create_boundary_particle_3D.F90"
#include "colloid_noslip.F90"
#include "colloid_compute_rotation_vector.F90"
#include "colloid_compute_rotation_matrix.F90"
#include "colloid_init_accumulation_matrix.F90"
#include "colloid_compute_accumulation_matrix.F90"
#include "colloid_compute_accumulation_vector.F90"
#include "colloid_collect_particles_interaction.F90"
#include "colloid_compute_interaction.F90"
#include "colloid_compute_lubrication_cc.F90"
#include "colloid_compute_lubrication_cw.F90"
#include "colloid_compute_repulsion_cc.F90"
#include "colloid_compute_repulsion_cw.F90"
#include "colloid_compute_interaction_implicit_all.F90"
#include "colloid_compute_interaction_implicit_velocity_pair.F90"
#include "colloid_compute_acceleration.F90"
#include "colloid_apply_body_force.F90"
#include "colloid_integrate_position.F90"
#include "colloid_integrate_velocity.F90"
#include "colloid_adjust_colloid.F90"
#include "colloid_compute_statistic.F90"
#include "colloid_polar_ellipse.F90"
#include "colloid_polar_star.F90"
#include "colloid_cartesian_ellipse.F90"
#include "colloid_spherical_ellipsoid.F90"
#include "colloid_cartesian_ellipsoid.F90"
#include "colloid_initialize_image.F90"
#include "colloid_compute_image.F90"
#include "colloid_nearest_image.F90"
#include "colloid_in_nearest_image.F90"
#include "colloid_in_relative_position.F90"
#include "colloid_particle_velocity.F90"
#include "colloid_find_force_extreme.F90"
#include "colloid_compute_dt_f.F90"
#include "colloid_set_flow_developed.F90"

      END MODULE Class_Colloid
      
      
