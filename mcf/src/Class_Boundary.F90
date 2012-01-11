      MODULE Class_Boundary
        !----------------------------------------------------
        ! Class	      : Boundary
        !----------------------------------------------------
        !
        !  Purpose    :
        !> \brief        Variables and corresponding operations
        !>               for boundary quantities.
        !>	   	
        !>               The variable memebers are private.
        !
        !  References   :
        !
        !  Remarks      :  
        !
        !
        !  Revisions    :  V0.1 29.09.2009, original version
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

        
        IMPLICIT NONE
        SAVE
        
        TYPE Boundary
           
           PRIVATE
           
           !-------------------------------------------------
           ! Member parameters
           !
           ! num_dim   : number of dimension.
           ! bcdef     : boundary definition.
           ! shear_rate: shear rate of different sides.
           ! shear_length: shear lenght for Lees-Edwards boundaires.
           ! shear_type : if there is wall, what type is it,
           !             1 : usuall wall with U0.
           !             2 : oscillating wall with U0*cos(w*t).
           ! shear_v0   : initial velocity for each side of the box.
           ! shear_v    : time-denpendent velocity for each side of the box.
           ! shear_freq : if oscillating shear, the frequency.
           ! rho_type    : 0 constant 1 calculated;
           ! noslip_type : no-slip conditons for wall
           !               1: Frozen particles.
           !               2: Morris J. P. et al. 1997
           ! dout      : minimal distance of a fluid particle 
           !             from wall, in case boundary particles
           !             are placed on the surface.
           ! drag      : exerted drag/force on the wall.
           ! drag_p    : pressure force.
           ! drag_v    : viscous force.
           ! drag_r    : random force.
           ! min_phys  : minimum boundaries of the domain.
           ! max_phys  : maximum boundaries of the domain.
           ! num_peri  : number of periodic boundaries.
           ! num_sym   : number of symmetry boundaires.
           ! num_wall_sym  : number of wall boundaires,
           !             handeled by PPM using symmetry.
           ! num_wall_solid : number of wall boundaires,
           !             handeled by MCF using solid particles.
           ! num_osci  : number of oscillating wall.
           ! num_le    : number of Lees-Edwards boundaries.
           !
           ! num_part_wall_solid : 
           !             number of wall boundary particles,
           !             created by MCF using solid particles.
           !-------------------------------------------------
           
           INTEGER                              :: num_dim           
           INTEGER, DIMENSION(:), POINTER       :: bcdef
           REAL(MK), DIMENSION(:,:), POINTER    :: shear_v0
           REAL(MK), DIMENSION(:,:), POINTER    :: shear_v
           REAL(MK), DIMENSION(:,:), POINTER    :: shear_rate
           REAL(MK), DIMENSION(:,:), POINTER    :: shear_length
           INTEGER, DIMENSION(:), POINTER       :: shear_type
           REAL(MK), DIMENSION(:), POINTER      :: shear_freq
           
           INTEGER                              :: rho_type
           INTEGER                              :: noslip_type
           REAL(MK)                             :: dout
           REAL(MK), DIMENSION(:,:),POINTER     :: drag
#ifdef __FORCE_SEPARATE
           REAL(MK), DIMENSION(:,:),POINTER     :: drag_p
           REAL(MK), DIMENSION(:,:),POINTER     :: drag_v
           REAL(MK), DIMENSION(:,:),POINTER     :: drag_r
#endif
           REAL(MK), DIMENSION(3)               :: min_phys
           REAL(MK), DIMENSION(3)               :: max_phys
           REAL(MK), DIMENSION(3)               :: min_phys_t
           REAL(MK), DIMENSION(3)               :: max_phys_t
         
           INTEGER                              :: num_peri
           INTEGER                              :: num_sym
           INTEGER                              :: num_wall_sym
           INTEGER                              :: num_wall_solid
           INTEGER                              :: num_osci
           INTEGER                              :: num_le
           INTEGER                              :: num_shear


           INTEGER                              :: num_part_wall_solid
           
           
        END TYPE Boundary
        
        INTERFACE boundary_new
           MODULE PROCEDURE boundary_init_default
           MODULE PROCEDURE boundary_init
        END INTERFACE

        INTERFACE boundary_noslip
           MODULE PROCEDURE boundary_noslip_solid
           MODULE PROCEDURE boundary_noslip_mirror
        END INTERFACE
        
      CONTAINS
        
#include "boundary_new.F90"
#include "boundary_finalize.F90"
#include "boundary_get.F90"
#include "boundary_set.F90"
#include "boundary_check_wall_solid_particle.F90"
#include "boundary_noslip_solid.F90"
#include "boundary_noslip_mirror.F90"
#include "boundary_collect_particles_interaction.F90"
#include "boundary_collect_colloid_interaction.F90"
#include "boundary_update_boundary.F90"
      
      END MODULE Class_Boundary
      
