      MODULE Class_Particles
        !----------------------------------------------------
        ! Class	      :	Particles
        !----------------------------------------------------
        !
        !  Purpose    :
        !> \brief       Variables and corresponding operations
        !>              for Particles infomation.
        !>	   	
        !>              The variable memebers are public 
        !  Remarks    :
        !
        !  References :
        !
        !  Revisions  : V0.1 03.03.2009 
        !
        !----------------------------------------------------
        !  Author     : Xin Bian
        !  Contact    : xin.bian@aer.mw.tum.de
        !  Dr. Marco Ellero's Emmy Noether Group,
        !  Prof. Dr. N. Adams' Chair of Aerodynamics,
        !  Faculty of Mechanical Engineering,
        !  Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        USE mcf_header
        USE Class_Debug
        USE Class_Tool
        
        USE Class_Control
        USE Class_Physics
        USE Class_Colloid
        USE Class_Rhs
        USE Class_StateEquation
        USE Class_Technique
        USE Class_Kernel
	USE Class_Random
        
        
        IMPLICIT NONE
        SAVE
        
        TYPE Particles
           
           TYPE(Control), POINTER               :: ctrl
           TYPE(Physics), POINTER               :: phys
           TYPE(Rhs), POINTER                   :: rhs
           TYPE(StateEquation), POINTER         :: stateEquation
           TYPE(Kernel), POINTER                :: kern
           TYPE(Technique), POINTER             :: tech
           LOGICAL                              :: pp_interact_cc
           LOGICAL                              :: pp_interact_cw
           
           !-------------------------------------------------
           ! General remarks:
           ! 
           ! Some variables/quantities are only useful when 
           ! some specific compiler flags are set. 
           ! However we define them all here,  and 
           ! initialize them in particles_new() and 
           ! finalize in particles_finalize().
           ! Only when they need to be given memory,
           ! we decide wethear they are used or allocated.
           ! This ways saves a lot of setting for compiler flags.
           !
           !
           ! x     : position
           ! v     : velocity
           ! rho   : mass / number density
           ! m     : mass
           ! p     : pressure
           ! id    : (1,:) particle ID
           !         (2,:) species  ID, i.e.,
           !          0 for fluid, > 0 for colloids
           !          < 0 for walls
           ! f     : force per unit mass / acceleration
           ! fa_min: minimum force acceleration.
           ! fa_max: maximum force acceleration.
           ! dt_f  : time step limit due to force acceleration.
           !
           ! fp    : pressure force per unit mass / acceleration
           ! fv    : viscous force per unit mass / acceleration
           ! fr    : randome force per unit mass / acceleration
           !           
           ! s     : stress tensor
           !
           !         s is calculated as rij*fij in the simulation,
           !         therefore, the real stress tensor with
           !         correct unit should be (m*rij*fij)/V,
           !         where V is volume of total fluid.
           !
           ! sp    : stress tensor from potential/pressure force
           ! sv    : stress tensor from viscous force
           ! sr    : stress tensor from random force.
           !
           ! Remark :
           !         Since PPM doesn't support tensor/matrix 
           !         for particles, i.e., it doesn't support 
           !         rank=3 dimensional data,
           !         we have to use dim**2 array for each particle, 
           !         if a matrix is needed.
           !
           !         In context of Fortran language, column 
           !         dominating convention is used, 
           !         i.e., index change first by rows then columns.
           !
           !
           ! s, sp, sv, sr
           ! matrix notation, 2D(3D):
           !  sxx   sxy  (sxz)
           !  syx   syy  (syz)
           ! (szx) (szy) (szz)
           !
           ! array notation of 2D in order:
           ! sxx, syx, sxy, syy;
           ! array notation of 3D in order:
           ! sxx, syx, szx, sxy, syy, szy,
           ! sxz, syz, szz.
           !
           ! u     : thermal energy.
           ! au    : acceleration of thermal energy.
           !
           ! Remark :
           !         1) Using PPM library, if there are
           !         ghosts particles, they are always
           !         in the end of the data struture
           !          continously.
           !-------------------------------------------------
           
           INTEGER                              :: num_dim
           REAL(MK)                             :: h
           REAL(MK), DIMENSION(:,:), POINTER    :: x
           REAL(MK), DIMENSION(:,:), POINTER    :: v
           REAL(MK), DIMENSION(:), POINTER      :: rho
           REAL(MK), DIMENSION(:), POINTER      :: rho_norm
           REAL(MK)                             :: rho_min
           REAL(MK)                             :: rho_max
           REAL(MK), DIMENSION(:), POINTER      :: m
           REAL(MK), DIMENSION(:), POINTER      :: p
           INTEGER, DIMENSION(:,:), POINTER     :: id
           REAL(MK), DIMENSION(:,:), POINTER    :: f

           REAL(MK)                             :: fa_min
           REAL(MK)                             :: fa_max
           REAL(MK)                             :: dt_f
           

           REAL(MK), DIMENSION(:,:), POINTER    :: fp
           REAL(MK), DIMENSION(:,:), POINTER    :: fv
           REAL(MK), DIMENSION(:,:), POINTER    :: fr           

           
           REAL(MK), DIMENSION(:,:), POINTER    :: s

           REAL(MK), DIMENSION(:,:), POINTER    :: sp
           REAL(MK), DIMENSION(:,:), POINTER    :: sv
           REAL(MK), DIMENSION(:,:), POINTER    :: sr
           
           
           REAL(MK), DIMENSION(:), POINTER      :: u
           REAL(MK), DIMENSION(:), POINTER      :: au
           
           !-------------------------------------------------
           ! Quantites more related to non-Newtonian fluid.
           !
           ! vgt   : velocity gradient tensor
           !
           ! evgt  : egenvector dynamics vgt
           ! eval  : egenvalues
           ! aeval : acceleration of egenvalues
           ! evec  : egenvectors
           ! aevec : acceleration of egenvectors           
           ! 
           ! ct    : conformation tensor
           ! act   : acceleation of conformation tensor
           !
           ! Remark :
           !         Since PPM doesn't support tensor/matrix 
           !         for particles, i.e., it doesn't support 
           !         rank=3 dimensional data,
           !         we have to use dim**2 array for each particle, 
           !         if a matrix is needed.
           !
           !         In context of Fortran language, column 
           !         dominating convention is used, 
           !         i.e., index change first by rows then columns.
           !
           ! vgt, evgt, ct or act's
           ! matrix notation, 2D(3D):
           !  Cxx   Cxy  (Cxz)
           !  Cyx   Cyy  (Cyz)
           ! (Czx) (Czy) (Czz)
           !
           ! array notation of 2D in order:
           ! Cxx, Cyx, Cxy, Cyy;
           ! array notation of 3D in order:
           ! Cxx, Cyx, Czx, Cxy, Cyy, Czy,
           ! Cxz, Cyz, Czz.
           !
           ! egenvalue 2D(3D):
           ! e1,e2(,e3).
           !
           ! egenvector matrix notation, 2D(3D):
           !
           !  ev1_x   ev2_x  (ev3_x)
           !  ev1_y   ev2_y  (ev3_y)
           ! (ev1_z) (ev2_z) (ev3_z)
           !
           ! 2D array notation in order :
           ! ev1_x, ev1_y, ev2_x, ev2_y
           ! 3D array notation in order :
           ! ev1_x, ev1_y, ev1_z, ev2_x, ev2_y, ev2_z,
           ! ev3_x, ev3_y, ev3_z
           !
           ! However, pressure tensor(pt) is not needed
           ! to communicate between processors and
           ! not using PPM,
           ! therefore is is matrix notation.
           !-------------------------------------------------
           
           REAL(MK), DIMENSION(:,:), POINTER    :: vgt
           REAL(MK), DIMENSION(:,:), POINTER    :: evgt
           REAL(MK), DIMENSION(:,:), POINTER    :: eval
           REAL(MK), DIMENSION(:,:), POINTER    :: aeval
           REAL(MK), DIMENSION(:,:), POINTER    :: evec
           REAL(MK), DIMENSION(:,:), POINTER    :: aevec
           REAL(MK), DIMENSION(:,:), POINTER    :: ct
           REAL(MK), DIMENSION(:,:), POINTER    :: act
           REAL(MK),DIMENSION(:,:,:),POINTER    :: pt
           
           ! Num of ids for each particle
           INTEGER                              :: num_id
           ! Particle id index
           INTEGER                              :: pid_idx
           ! Spieces id index
           INTEGER                              :: sid_idx
           
           !----------------------------------------------------
           ! real :  real particles on local process;
           ! all  :  all particles including ghost;
           ! ghost:  ghost particles.
           !
           ! all = real + ghost
           !
           ! fluid   : number of real fluid particles;
           ! sym     : number of symmetry boundary particles.
           !           ( are all ghosts particles)
           ! wall_sym: number of wall boundary particles created
           !           by PPM using symmetry/mirror particles.
           !           ( are all ghosts particles)
           ! wall_solid: inital number of wall boundary real particles created
           !           by MCF using solid particles.
           ! wall_solid_real : number of wall boundary ghost particles
           !           created by MCF using solid particles.
           !           ( are real particles)
           ! wall_solid_ghost : number of wall boundary ghost particles
           !           created by MCF using solid particles.
           !           ( are ghosts particles)
           ! le      : number of Lees-Edwards boundary particles.
           !           ( are ghost particles)
           !
           ! shear   : number of wall boundary particles,
           !           created by ppm using symmetry/mirror,
           !           or by MCF as solid wall, or number of
           !           Lees-Edwards boundary particles.
           !
           !
           ! colloid : number of colloid boundary particles.
           !
           ! real = fluid + wall_solid + colloid
           !-----------------------------------------------------
           
           INTEGER                              :: num_part_real
           INTEGER                              :: num_part_all
           INTEGER                              :: num_part_ghost
           
           INTEGER                              :: num_part_fluid
           INTEGER                              :: num_part_sym
           INTEGER                              :: num_part_wall_sym
           INTEGER                              :: num_part_wall_solid
           INTEGER                              :: num_part_wall_solid_real
           INTEGER                              :: num_part_wall_solid_ghost
           INTEGER                              :: num_part_le
           INTEGER                              :: num_part_colloid
           
           !----------------------------------------------------
           ! Record the index and species ID of each :
           !
           ! symmetry boundary particle;
           ! ( are all ghost particles)
           ! wall_using symmetry boundary particles;
           ! ( are all ghost particles)
           ! wall using solid boundary particles;
           ! ( part of are real, part of are ghost particles)
           ! Lees-Edwards boundary particles.
           ! ( are all ghost particles)
           ! colloid boundary particle;
           ! ( part of are real, part of are ghost particles)
           !----------------------------------------------------
           
           INTEGER, DIMENSION(:,:), POINTER     :: part_sym_list
           INTEGER, DIMENSION(:,:), POINTER     :: part_walL_sym_list
           INTEGER, DIMENSION(:,:), POINTER     :: part_wall_solid_real_list
           INTEGER, DIMENSION(:,:), POINTER     :: part_wall_solid_ghost_list
           INTEGER, DIMENSION(:,:), POINTER     :: part_le_list
           INTEGER, DIMENSION(:,:), POINTER     :: part_colloid_list
           
           TYPE(Tool)                           :: tool
           TYPE(Random), POINTER                :: random
           
        END TYPE Particles
        
        
        INTERFACE particles_new
           MODULE PROCEDURE particles_init
        END INTERFACE
        
      CONTAINS
        
#include "particles_new.F90"
#include "particles_finalize.F90"
#include "particles_get.F90"
#include "particles_set.F90"
#include "particles_reset.F90"
#include "particles_init_global_inter.F90"
#include "particles_init_global_exter.F90"
#include "particles_init_global_assign_id.F90"
#include "particles_init_partial_inter.F90"
#include "particles_init_partial_exter.F90"
#include "particles_set_colloid_on_lattice.F90"
#include "particles_decompose_global.F90"
#include "particles_decompose_partial.F90"
#include "particles_map_ghost_get.F90"
#include "particles_map_ghost_put.F90"
#include "particles_compute_mass.F90"
#include "particles_compute_density.F90"
#include "particles_normalize_density.F90"
#include "particles_find_density_extreme.F90"
#include "particles_compute_pressure.F90"
#include "particles_compute_interaction.F90"
#include "particles_find_force_extreme.F90"
#include "particles_compute_dt_f.F90"
#include "particles_apply_body_force.F90"
#include "particles_compute_act.F90"
#include "particles_compute_evgt.F90"
#include "particles_compute_aeval.F90"
#include "particles_compute_aevec.F90"
#include "particles_compute_ct.F90"
#include "particles_compute_pressure_tensor.F90"
#include "particles_collect_colloid_interaction.F90"
#include "particles_collect_boundary_interaction.F90"
#include "particles_integrate_position.F90"
#include "particles_integrate_velocity.F90"
#include "particles_compute_colloid_relative_position.F90"
#include "particles_compute_colloid_absolute_position.F90"
#include "particles_integrate_boundary_position.F90"
#include "particles_adjust_particles.F90"
#include "particles_integrate_eval.F90"
#include "particles_integrate_evec.F90"
#include "particles_integrate_ct.F90"
#include "particles_integrate_potential_energy.F90"
#include "particles_set_boundary_ghost_id.F90"
#include "particles_set_boundary_velocity.F90"
#include "particles_set_boundary_ghost_velocity.F90"
#include "particles_reset_boundary_velocity.F90"
#include "particles_reset_boundary_interaction.F90"
#include "particles_reset_boundary_ghost_interaction.F90"
#include "particles_set_colloid_velocity.F90"
#include "particles_reset_colloid_velocity.F90"
#include "particles_reset_colloid_interaction.F90"
#include "particles_set_flow_developed.F90"

      END MODULE Class_Particles
      
      
