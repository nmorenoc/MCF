      SUBROUTINE colloid_noslip_Morris_ellipse(this, &
           xf,xc,vf,vc,sid_c,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_noslip_Morris_ellipse
        !----------------------------------------------------
        ! Purpose     : For 2D ellipse.
        !
        ! Revision    : V0.1 16.03.2009, orignal version.
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
        ! Local variables start here :
        ! 
        ! i     : index.
        ! dim   : number of dimension.
        ! xcoll : position of the center of the colloid.
        ! vcoll : velocity of the center of the colloid.
        ! r_xf    : relative position from fluid particle 
        !           to the center.
        ! r_xc    : relative position from colloid boundary
        !           particle to center.
        ! r_vc    : relative velocity of colloid boundary 
        !           particle to the center.
        ! theta   : angel of fluid particle to x+ direction.
        ! d_scoll : distance from surface to the center.     
        ! d_fcoll : distance from fluid particle to the center.
        ! d_ccoll : distance from the colloid boundary particle
        !           to the center.
        ! d_fs    : distance from fluid to tangent surface,
        !           i.e., shortes distance.
        ! xs      : position on the surface which has
        !           shortestd distance to fluid particle.
        !           (This is the relative position to
        !            the center.)
        ! nvector : unit vector pointing from xs to
        !           r_xf.
        ! d_cs    : distance from xc to tangent surface.
        ! corr    : d_cs/d_fs.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: i
        INTEGER                                 :: dim
        REAL(MK), DIMENSION(3)                  :: xcoll
        REAL(MK), DIMENSION(3)                  :: vcoll
        REAL(MK), DIMENSION(3)                  :: r_xf
        REAL(MK), DIMENSION(3)                  :: r_xc
        REAL(MK), DIMENSION(3)                  :: r_vc
        REAL(MK)                                :: theta
        REAL(MK)                                :: d_scoll
        REAL(MK)                                :: d_fcoll
        REAL(MK)                                :: d_ccoll
        REAL(MK), DIMENSION(3)                  :: xs
        REAL(MK)                                :: d_fs
        REAL(MK), DIMENSION(3)                  :: nvector
        REAL(MK)                                :: d_cs
        REAL(MK)                                :: corr
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        dim = this%num_dim
        r_xc(:) = 0.0_MK
        r_vc(:) = 0.0_MK     
        
     
        !----------------------------------------------------
        ! Get the nearest image colloid center's to the
        ! boundary particle.
        !----------------------------------------------------
        
        CALL colloid_in_nearest_image(this,xc(1:dim),sid_c, &
             xcoll(1:dim),r_xc(1:dim),vcoll(1:dim),stat_info_sub)
      
        !----------------------------------------------------
        ! Get relative position of fluid particle to the center.
        !----------------------------------------------------
        
        r_xf(1:dim) = xf(1:dim) - xcoll(1:dim)
        
        !----------------------------------------------------
        ! Get translational angular velocity if rotating.
        !----------------------------------------------------
        
        r_vc(:) = 0.0_MK
        
        IF ( this%rotate ) THEN
           
           CALL tool_cross_product(this%tool,&
                this%omega(1:3,sid_c), r_xc(1:3),&
                r_vc(1:3),stat_info_sub)
           
        END IF

        !----------------------------------------------------
        ! Angle between r_xf and x+ direction. Get at this 
        ! angle what is the distance from the surface.
        !----------------------------------------------------
        
        theta = polar_angle(r_xf(1),r_xf(2))
        
        d_scoll = polar_ellipse_r(this%radius(1,sid_c), &
             this%radius(2,sid_c),theta,this%theta(3,sid_c))
        
        !----------------------------------------------------
        ! Distance between a fluid particle and the center.
        !----------------------------------------------------
        
        d_fcoll = 0.0_MK
        DO i = 1, dim
           d_fcoll = d_fcoll + r_xf(i)**2
        END DO
        d_fcoll = SQRT(d_fcoll)
        
        !----------------------------------------------------
        ! Check if penetration happens, i.e., a fluid 
        ! particle goes inside a colloid object. If yes,
        ! simply assign the same velocity as the colloid.
        !----------------------------------------------------
        
        IF( d_fcoll <= d_scoll ) THEN
           
           vc(1:dim) = vcoll(1:dim) + r_vc(1:dim)
           
           !PRINT *, "Fluid particle inside the colloid sphere!"          
           !stat_info = -1
           GOTO 9999
           
        END IF
        
        !----------------------------------------------------
        ! Angle between r_xc and x+ direction. Get at this 
        ! angle what is the distance from the surface.
        !----------------------------------------------------
        
        theta = polar_angle(r_xc(1),r_xc(2))
        
        d_scoll = polar_ellipse_r(this%radius(1,sid_c), &
             this%radius(2,sid_c),theta,this%theta(3,sid_c) )
        
        !----------------------------------------------------
        ! Distance between colloid boundary particle 
        ! and the center of a colloid.
        !----------------------------------------------------
        
        d_ccoll = SQRT(DOT_PRODUCT(r_xc(1:dim),r_xc(1:dim)))
        
        !----------------------------------------------------
        ! Inconsistent movement happens, if colloid boundary 
        ! particle goes out of the geometry.
        !----------------------------------------------------
        
        IF( d_ccoll > d_scoll + this%dout ) THEN
           PRINT *, "colloid_noslip_Morris_ellipse : ", &
                "Colloid particle goes out of colloid ellipse !"
           PRINT *, "dist,ra,rb,xc,xcolloid,xcoll :", &
                d_ccoll,this%radius(1,sid_c),this%radius(2,sid_c), &
                xc(1:dim),this%x(1:dim,sid_c),xcoll(1:dim)
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Find the shortest distance from fluid particle
        ! to the surface of ellipse.
        ! s_x(:) contains the surface point.
        !----------------------------------------------------
        
        CALL cartesian_ellipse_shortestD(this%radius(1,sid_c), &
             this%radius(2,sid_c),this%theta(3,sid_c), &
             r_xf(1),r_xf(2),xs(1),xs(2),d_fs)
        
        !----------------------------------------------------
        ! Normal vector pointing from xs to r_xf.
        !----------------------------------------------------
        
        nvector(1:dim) = &
             ( r_xf(1:dim)- xs(1:dim) ) /  d_fs
        
        !----------------------------------------------------
        ! Map (xs-r_xc) onto nvector and calculate length.
        !----------------------------------------------------
        
        d_cs = 0.0_MK
        DO i =1, dim           
           d_cs = d_cs+ (xs(i)-r_xc(i)) * nvector(i)
        END DO
        
        !----------------------------------------------------
        ! If the fuild particle lies exactly on
        ! the surface, assign a minimum distance.
        !----------------------------------------------------
        
        IF( d_fs < ABS( mcf_machine_zero) ) THEN
           
           d_fs = mcf_machine_zero
           
        END IF
        
        !----------------------------------------------------
        ! Ratio of the distances, used for extrapolation.
        !----------------------------------------------------
        
        corr = d_cs / d_fs
        
        !----------------------------------------------------
        ! Set the maximum ratio of the extrapolation.
        !----------------------------------------------------
        
        IF (corr > 0.5_MK) THEN
           
           corr = 0.5_MK
           
        END IF
        
        !----------------------------------------------------
        ! Extrapolate velocity for the colloid boundary 
        ! particle, considering the movment of colloid.
        !----------------------------------------------------
        
        vc(1:dim) = -corr *&
             (vf(1:dim)-vcoll(1:dim)-r_vc(1:dim)) + &
             vcoll(1:dim) + r_vc(1:dim)
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_noslip_Morris_ellipse
