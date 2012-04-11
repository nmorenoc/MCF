      SUBROUTINE colloid_noslip_Morris(this,xf,xc,vf,vc,sid_c,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_noslip_Morris
        !----------------------------------------------------
        !
        ! Purpose     : Calculate the extrapolated velocity
        !               for a colloid boundary particle
        !               in order to have no slip condition
        !               on the surface of a colloid.
        !
        ! Remark      : Extend no slip condition from 
        !               Morris J.P. et al. 1997.
        !               Check Bian et. al 2012 for details.
        !
        ! Revision    : V0.3 23.03.2010,
        !               including star-like shape.
        !
        !               V0.2 27.11.2009, including
        !               Lees-Edwards boundary.
        !
        !               V0.1 01.03.2009, orignal version.
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
        ! this  : object of a colloid.
        ! xf    : position of a fluid particle.
        ! xc    : position of a colloid boundary particle.
        ! vf    : velocity of a fluid particle.
        ! sid_c : species ID of a colloid boundary particle.
        !
        ! Output
        !
        ! vc    : extrapolated velocity for the colloid
        !         boundary particle.
        ! stat_info : status of the routine.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xf
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xc
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vf
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: vc
        INTEGER, INTENT(IN)                     :: sid_c        
        INTEGER, INTENT(OUT)                    :: stat_info 
        
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        !----------------------------------------------------
        ! Check the shape of a colloid and then choose
        ! sub-routine correspondingly.
        !----------------------------------------------------
        
        IF ( this%num_dim == 2 ) THEN

           SELECT CASE( this%shape(sid_c) )
              
           CASE (mcf_colloid_shape_cylinder)
              
              !----------------------------------------------
              ! For a 2D cylinder
              !----------------------------------------------
              
              CALL colloid_noslip_Morris_cylinder_2D(this, &
                   xf,xc,vf,vc,sid_c,stat_info_sub)
              
              IF( stat_info_sub /= 0 ) THEN
                 PRINT *, __FILE__, ":", __LINE__
                 stat_info = -1
                 GOTO 9999
              END IF
              
           CASE (mcf_colloid_shape_disk)
              
              !----------------------------------------------
              ! For a 2D disk
              !----------------------------------------------
              
              CALL colloid_noslip_Morris_disk(this, &
                   xf,xc,vf,vc,sid_c,stat_info_sub)
              
              IF( stat_info_sub /= 0 ) THEN
                 PRINT *, __FILE__, ":", __LINE__
                 stat_info = -1
                 GOTO 9999
              END IF
              
           CASE (mcf_colloid_shape_ellipse)
              
              !----------------------------------------------
              ! For a 2D ellipse.
              !----------------------------------------------
              
              CALL colloid_noslip_Morris_ellipse(this, &
                   xf,xc,vf,vc,sid_c,stat_info_sub)
              
              IF( stat_info_sub /= 0 ) THEN
                 PRINT *, __FILE__, ":", __LINE__
                 stat_info = -1
                 GOTO 9999
              END IF
#if 0              
           CASE (mcf_colloid_shape_dicolloid)
              
              !----------------------------------------------
              ! For a 2D dicolloid
              !----------------------------------------------
              
              CALL colloid_noslip_Morris_dicolloid_2D(this, &
                   xf,xc,vf,vc,sid_c,stat_info_sub)
              
              IF( stat_info_sub /= 0 ) THEN
                 PRINT *, __FILE__, ":", __LINE__
                 stat_info = -1
                 GOTO 9999
              END IF
#endif         
           CASE (mcf_colloid_shape_star)
              
              !----------------------------------------------
              ! For a 2D star-like shape (penut, star...)
              !----------------------------------------------
              
              CALL colloid_noslip_Morris_star_2D(this, &
                   xf,xc,vf,vc,sid_c,stat_info_sub)
              
              IF( stat_info_sub /= 0 ) THEN
                 PRINT *, __FILE__, ":", __LINE__
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END SELECT ! shape
           
        ELSE IF ( this%num_dim == 3 ) THEN
           
           SELECT CASE( this%shape(sid_c) )
              
           CASE (mcf_colloid_shape_cylinder)
              
              !----------------------------------------------
              ! For a 3D cylinder
              !----------------------------------------------
              
              CALL colloid_noslip_Morris_cylinder_3D(this, &
                   xf,xc,vf,vc,sid_c,stat_info_sub)
              
              IF( stat_info_sub /= 0 ) THEN
                 PRINT *, __FILE__, ":", __LINE__
                 stat_info = -1
                 GOTO 9999
              END IF

           CASE (mcf_colloid_shape_sphere)
              
              !----------------------------------------------
              ! For a 3D sphere.
              !----------------------------------------------
              
              CALL colloid_noslip_Morris_sphere(this, &
                   xf,xc,vf,vc,sid_c,stat_info_sub)
              
               IF( stat_info_sub /= 0 ) THEN
                 PRINT *, __FILE__, ":", __LINE__
                 stat_info = -1
                 GOTO 9999
              END IF
              
           CASE (mcf_colloid_shape_ellipsoid)
              
              !----------------------------------------------
              ! For a 3D ellipsoid.
              !----------------------------------------------
              
              CALL colloid_noslip_Morris_ellipsoid(this, &
                   xf,xc,vf,vc,sid_c,stat_info_sub)
              
              IF( stat_info_sub /= 0 ) THEN
                 PRINT *, __FILE__, ":", __LINE__
                 stat_info = -1
                 GOTO 9999
              END IF
             
           CASE (mcf_colloid_shape_dicolloid)
              
              !----------------------------------------------
              ! For a 3D dicolloid
              !----------------------------------------------
              
              CALL colloid_noslip_Morris_dicolloid_3D(this, &
                   xf,xc,vf,vc,sid_c,stat_info_sub)
              
              IF( stat_info_sub /= 0 ) THEN
                 PRINT *, __FILE__, ":", __LINE__
                 stat_info = -1
                 GOTO 9999
              END IF
        
           END SELECT ! shape
           
        END IF
           
9999       CONTINUE
           
        RETURN
        
      END SUBROUTINE colloid_noslip_Morris

#include "colloid_noslip_Morris_cylinder_2D.F90"
#include "colloid_noslip_Morris_cylinder_3D.F90"
#include "colloid_noslip_Morris_disk.F90"
#include "colloid_noslip_Morris_sphere.F90"
#include "colloid_noslip_Morris_ellipse.F90"
#include "colloid_noslip_Morris_ellipsoid.F90"
#include "colloid_noslip_Morris_dicolloid_3D.F90"
#include "colloid_noslip_Morris_star_2D.F90"
