      SUBROUTINE colloid_noslip_Morris_ellipse(this, &
           xf,xc,vf,vc,sid_c,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_noslip_Morris_ellipse
        !----------------------------------------------------
        ! Purpose     : For 2D ellipse.
        !
        ! Revision    : V0.2 14.03.2012, make sure that
        !               the extrapolated velocity of the
        !               boundary particle makes the surface
        !               point velocity with noslip condition.
        !               
        !               V0.1 16.03.2009, orignal version.
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
        ! vc        : extrapolated velocity for the colloid
        !             boundary particle.
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
        ! Local variables start here :
        ! 
        ! i       : index.
        ! dim     : number of dimension.
        ! xcoll   : position of the center of the colloid.
        ! vcoll   : velocity of the center of the colloid.
        ! r_xf    : relative position from fluid particle 
        !           to the center.
        ! r_xc    : relative position from colloid boundary
        !           particle to center.
        ! r_vc    : relative velocity of colloid boundary 
        !           particle to the center.
        ! theta   : angle of fluid particle to x+ direction.
        ! d_scoll : distance from surface to the center.     
        ! d_fcoll : distance from fluid particle to the center.
        ! d_ccoll : distance from the colloid boundary particle
        !           to the center.
        ! d_fs    : distance from fluid to tangent surface,
        !           i.e., shortest distance.
        ! r_xs     : position on the surface which has
        !           shortest distance to fluid particle.
        !           (This is the relative position to
        !           the center.)
        ! nvector : unit vector pointing from r_xs to
        !           r_xf.
        ! d_cs    : distance from r_xc to tangent surface.
        ! corr    : d_cs/d_fs.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        REAL(MK), DIMENSION(3)                  :: xcoll
        REAL(MK), DIMENSION(3)                  :: vcoll
        REAL(MK), DIMENSION(3)                  :: r_xf
        REAL(MK), DIMENSION(3)                  :: r_xc
        REAL(MK), DIMENSION(3)                  :: r_vc
        REAL(MK), DIMENSION(3)                  :: r_xs
        REAL(MK), DIMENSION(3)                  :: r_vs
#if 0
        REAL(MK)                                :: theta
        REAL(MK)                                :: d_scoll
        REAL(MK)                                :: d_fcoll
        REAL(MK)                                :: d_ccoll
#endif
        REAL(MK)                                :: d_fs
        REAL(MK), DIMENSION(3)                  :: nvector
        REAL(MK)                                :: d_cs
        REAL(MK)                                :: corr
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        dim     = this%num_dim
        r_xc(:) = 0.0_MK
        r_vc(:) = 0.0_MK     
        
        !----------------------------------------------------
        ! Get the nearest image colloid center's to the
        ! boundary particle.
        !----------------------------------------------------
        
        CALL colloid_in_nearest_image(this,xc(1:dim),sid_c,&
             xcoll(1:dim),r_xc(1:dim),vcoll(1:dim),&
             stat_info_sub)
      
        !----------------------------------------------------
        ! Get relative position of fluid particle to the center.
        !----------------------------------------------------
        
        r_xf(1:dim) = xf(1:dim) - xcoll(1:dim)
        

#if 0
        !----------------------------------------------------
        ! Check if the fluid particle is inside of 
        ! the geometry.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Angle between r_xf and x+ direction. Get at this 
        ! angle, what is the distance from origin to 
        ! the surface.
        ! The distance d_scoll is used for penetration check.
        !----------------------------------------------------
        
        theta   = colloid_polar_angle(r_xf(1),r_xf(2))
        d_scoll = colloid_polar_ellipse_r(this%radius(1,sid_c), &
             this%radius(2,sid_c),theta,this%theta(3,sid_c))
        
        !----------------------------------------------------
        ! Distance between a fluid particle and the center.
        !----------------------------------------------------
        
        d_fcoll = SQRT(DOT_PRODUCT(r_xf(1:dim),r_xf(1:dim)))
        
        !----------------------------------------------------
        ! Check if penetration happens, i.e., a fluid 
        ! particle goes inside a colloid object. 
        ! If yes,
        ! simply assign the same velocity as the colloid.
        !----------------------------------------------------
        
        IF( d_fcoll <= d_scoll ) THEN
           
           !-------------------------------------------------
           ! Assign traslational velocity.
           !-------------------------------------------------
           
           vc(1:dim) = vcoll(1:dim)
           
           !-------------------------------------------------
           ! If rotating, add relative velocity of the 
           ! boundary particle to its center of geometry.
           ! i.e., equavilent translational 
           ! velcoity from angular velocity.
           !-------------------------------------------------
           
           IF ( this%rotate ) THEN
           
              r_vc(:) = 0.0_MK
              
              CALL tool_cross_product(this%tool,&
                   this%omega(1:3,sid_c,1), r_xc(1:3),&
                   r_vc(1:3),stat_info_sub)
           
              vc(1:dim) = vc(1:dim) + r_vc(1:dim)
              
           END IF
           
           !PRINT *, "Fluid particle inside the geometry!"
           !stat_info = -1
           GOTO 9999
           
        END IF
        
        !----------------------------------------------------
        ! Check if the boundary particle is outside of
        ! the geometry.
        !----------------------------------------------------
      
        !----------------------------------------------------
        ! Angle between r_xc and x+ direction. Get at this 
        ! angle what is the distance from origin to 
        ! the surface.
        ! The distance d_scoll is used for consistent moving
        ! check.        
        !----------------------------------------------------
        
        theta   = colloid_polar_angle(r_xc(1),r_xc(2))
        d_scoll = colloid_polar_ellipse_r(this%radius(1,sid_c), &
             this%radius(2,sid_c),theta,this%theta(3,sid_c) )
        
        !----------------------------------------------------
        ! Distance between colloid boundary particle 
        ! and the center of a colloid.
        !----------------------------------------------------
        
        d_ccoll = SQRT(DOT_PRODUCT(r_xc(1:dim),r_xc(1:dim)))
        
        !----------------------------------------------------
        ! Inconsistent movement happens, if the boundary 
        ! particle goes out of the geometry.
        !----------------------------------------------------
        
        IF( d_ccoll > d_scoll + this%dout ) THEN
           PRINT *, __FILE__, ":", __LINE__, &
                "boundary particle goes out of geometry !"
           PRINT *, "dist,xc,xcolloid,xcoll,theta :", &
                d_ccoll,xc(1:dim),&
                this%x(1:dim,sid_c),xcoll(1:dim), &
                this%theta(3,sid_c)
           
           stat_info = -1
           GOTO 9999
        END IF
        
#endif
        
        !----------------------------------------------------
        ! Find the shortest distance from the fluid particle
        ! to the surface of the ellipse.
        ! r_xs(:) contains the surface point.
        !----------------------------------------------------
        
        CALL colloid_cartesian_ellipse_shortestD(&
             this%radius(1,sid_c),this%radius(2,sid_c),&
             this%theta(3,sid_c),&
             r_xf(1),r_xf(2),r_xs(1),r_xs(2),&
             d_fs,stat_info_sub)
        
        IF ( stat_info_sub /=0 ) THEN
           PRINT *, __FILE__, ":", __LINE__
           stat_info = -1
           GOTO 9999
        END IF

        !----------------------------------------------------
        ! If the fuild particle lies exactly on
        ! the surface, assign the maximum ratio.
        !----------------------------------------------------
        
        IF( d_fs < ABS( mcf_machine_zero) ) THEN
           
           corr = mcf_colloid_dist_ratio_max
           
        ELSE
           
           !-------------------------------------------------
           ! Normal vector pointing from r_xs to r_xf.
           !-------------------------------------------------
           
           nvector(1:dim) = &
                ( r_xf(1:dim)- r_xs(1:dim) ) /  d_fs
           
           !-------------------------------------------------
           ! Map (r_xs-r_xc) onto nvector and 
           ! calculate it length.
           !-------------------------------------------------
        
           d_cs = DOT_PRODUCT((r_xs(1:dim)-r_xc(1:dim)), nvector(1:dim))
           
           
           !-------------------------------------------------
           ! Ratio of the distances, used for extrapolation.
           !-------------------------------------------------
           
           corr = d_cs / d_fs
           
           !-------------------------------------------------
           ! Set the maximum ratio of the extrapolation.
           !-------------------------------------------------
           
           IF ( corr > mcf_colloid_dist_ratio_max ) THEN
              
              corr = mcf_colloid_dist_ratio_max
              
           END IF
           
        END IF ! d_fs == 0
        
        !----------------------------------------------------
        ! Extrapolate velocity for the colloid boundary 
        ! particle, considering the movment of colloid.
        !----------------------------------------------------
        
        r_vs(:) = 0.0_MK
        
        CALL tool_cross_product(this%tool,&
             this%omega(1:3,sid_c,1), r_xs(1:3),&
             r_vs(1:3),stat_info_sub)
        
        vc(1:dim) = -corr *&
             (vf(1:dim)-vcoll(1:dim)-r_vs(1:dim)) + &
             vcoll(1:dim) + r_vs(1:dim)
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_noslip_Morris_ellipse
      
