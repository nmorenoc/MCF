      SUBROUTINE colloid_noslip_Morris_star(this, &
           xf,xc,vf,vc,sid_c,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_noslip_Morris_star
        !----------------------------------------------------
        ! Purpose     : For 2D star velocity no-slip
        !               boundary condition.
        !
        ! Reference   : Bian et al. 2011.
        !
        ! Revision    : V0.2 04.03.2011, new version according
        !               to Bian et al. 2011.
        !              
        !               V0.1 23.03.2010, orignal version.
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
        ! i       : index.
        ! dim     : number of dimension.
        !
        ! #nearest image of center of the colloid to xc.#
        ! xcoll   : position of the nearest image center
        !           of the colloid to particle c.
        ! vcoll   : velocity of the image center.
        !
        ! r_xf    : relative position from fluid particle 
        !           f to the iamge center.
        ! r_xc    : relative position from boundary particle
        !           c to the image center.
        ! r_xs    : relative position from surface point 
        !           s to the image center.
        ! r_vc    : relative velocity of boundary particle 
        !           c to the image center.
        ! r_vs    : relative velocity of surface point 
        !           s to the image center.
        ! d_fcoll : distance from f to the image center.
        ! d_scoll : distance from s to the image center.
        ! d_ccoll : distance from c to the image center.
        ! d_fs    : distance from f to tangent surface.
        ! d_u     : d_cs/d_fs+1.
        ! corr    : d_cs/d_fs.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: i
        INTEGER                                 :: dim
        REAL(MK), DIMENSION(3)                  :: xcoll
        REAL(MK), DIMENSION(3)                  :: vcoll
        REAL(MK)                                :: theta        
        REAL(MK), DIMENSION(3)                  :: r_xf
        REAL(MK), DIMENSION(3)                  :: r_xc
        REAL(MK), DIMENSION(3)                  :: r_xs
        REAL(MK), DIMENSION(3)                  :: r_vc
        REAL(MK), DIMENSION(3)                  :: r_vs
        REAL(MK)                                :: d_fcoll
        REAL(MK)                                :: d_scoll
        REAL(MK)                                :: d_ccoll
        REAL(MK)                                :: d_u
        REAL(MK)                                :: d_fs
        REAL(MK)                                :: corr
        
        REAL(MK), DIMENSION(3)                  :: r_xp
        LOGICAL                                 :: convex

        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        dim = this%num_dim
        
        xcoll(:)   = 0.0_MK
        vcoll(:)   = 0.0_MK
        r_xf(:)    = 0.0_MK
        r_xc(:)    = 0.0_MK
        r_xs(:)    = 0.0_MK
        r_vc(:)    = 0.0_MK
        r_vs(:)    = 0.0_MK        
        
        !----------------------------------------------------
        ! Get the nearest image of colloid center's to the
        ! boundary particle c.
        !----------------------------------------------------
        
        CALL colloid_in_nearest_image(this,xc(1:dim),sid_c,&
             xcoll(1:dim),r_xc(1:dim),vcoll(1:dim),stat_info_sub)
        
        
        !----------------------------------------------------
        ! Get relative position of fluid particle to 
        ! the image center.
        !----------------------------------------------------
        
        r_xf(1:2) = xf(1:2) - xcoll(1:2)
        
        !----------------------------------------------------
        ! Get equavilent translational velocity of the
        ! angular velocity, if it is rotating.
        !----------------------------------------------------
        
        IF ( this%rotate ) THEN
           
           CALL tool_cross_product(this%tool,&
                this%omega(1:3,sid_c), r_xc(1:3),&
                r_vc(1:3),stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Distance between f and the image center.
        !----------------------------------------------------
        
        d_fcoll = SQRT(DOT_PRODUCT(r_xf(1:2),r_xf(1:2)))
        
        !----------------------------------------------------
        ! Check if penetration happens, i.e., a fluid 
        ! particle goes inside a colloid object. If yes, 
        ! simply assign the same velocity as the colloid.
        !----------------------------------------------------
        
        theta = polar_angle(r_xf(1), r_xf(2))
        
        d_scoll = &
             polar_star_r(this%radius(1,sid_c),this%radius(2,sid_c), &
             REAL(this%freq(sid_c),MK),theta,this%theta(3,sid_c))
        
        IF( d_fcoll <= d_scoll ) THEN
           
           vc(1:2) = vcoll(1:2) + r_vc(1:2)
           !PRINT *, "Penetration of f is not neccessary an error!" 
           !PRINT *, "dist,radius,xf,xcolloid,xcoll :", &
           !     d_fcoll,this%radius(1,sid_c),xf(1:dim),&
           !     this%x(1:dim,sid_c),xcoll(1:dim)       
           !stat_info = -1
           GOTO 9999
           
        END IF
        
        !----------------------------------------------------
        ! Distance between c and the image center.
        !----------------------------------------------------
        
        d_ccoll = SQRT(DOT_PRODUCT(r_xc(1:2),r_xc(1:2)))
       
        !----------------------------------------------------
        ! Inconsistent movement happens, if c 
        ! goes out of the geometry.
        !----------------------------------------------------
        
        theta = polar_angle(r_xc(1), r_xc(2))        
        
        d_scoll = &
             polar_star_r(this%radius(1,sid_c),this%radius(2,sid_c), &
             REAL(this%freq(sid_c),MK),theta,this%theta(3,sid_c))
        
        IF( d_ccoll > d_scoll + this%dout ) THEN
           PRINT *, "colloid_noslip_Morris_star : ", &
                "Colloid particle goes out of colloid star !"
           PRINT *, "dist,ra,rb,xc,xcolloid,xcoll :", &
                d_ccoll,this%radius(1,sid_c),this%radius(2,sid_c), &
                xc(1:dim),this%x(1:dim,sid_c),xcoll(1:dim)
           stat_info = -1
           GOTO 9999
        END IF
        
        convex = .TRUE.
        
#if 0
        CALL polar_star_intersectP(this%radius(1,sid_c),&
             this%radius(2,sid_c),REAL(this%freq(sid_c),MK), &
             this%theta(3,sid_c),r_xc(1),r_xc(2), &
             r_xf(1),r_xf(2),r_xp(1),r_xp(2),stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "colloid_noslip_Morris_star : ", &
                "Finding intersecting point has problem"
           stat_info = -1
           GOTO 9999
        END IF
        
        convex = polar_star_convex(REAL(this%freq(sid_c),MK),&
             this%theta(3,sid_c),r_xp(1),r_xp(2),stat_info_sub)
        
#endif
        !----------------------------------------------------
        ! Get the minimal distance between f and the surface.
        !----------------------------------------------------
        
        IF ( convex ) THEN
           
           CALL polar_star_shortestD(this%radius(1,sid_c),&
                this%radius(2,sid_c),REAL(this%freq(sid_c),MK), &
                this%theta(3,sid_c),r_xf(1),r_xf(2), &
                r_xs(1),r_xs(2),d_fs,stat_info_sub)
           
        ELSE
           
           CALL polar_star_shortestD(this%radius(1,sid_c),&
                this%radius(2,sid_c),REAL(this%freq(sid_c),MK), &
                this%theta(3,sid_c),r_xc(1),r_xc(2), &
                r_xs(1),r_xs(2),d_fs,stat_info_sub)
           
        END IF
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "colloid_noslip_Morris_star : ", &
                "Finding shortest distance has problem"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! If the fuild particle lies exactly on
        ! the surface, i.e., f is same as s,
        ! assign a minimum ratio.
        !----------------------------------------------------
        
        IF ( ABS(r_xs(1)-r_xf(1)) < mcf_machine_zero .AND. &
             ABS(r_xs(2)-r_xf(2)) < mcf_machine_zero ) THEN
           
           corr = 0.5_MK
           
        ELSE
           
           !-------------------------------------------------
           ! Map vector cf onto sf(get nf), 
           ! then divided by sf^2,
           ! therefore d_u=(d_c+df)/df.
           !-------------------------------------------------
           
           d_u = ((r_xc(1)-r_xf(1))*(r_xs(1)-r_xf(1)) + &
                (r_xc(2)-r_xf(2))*(r_xs(2)-r_xf(2))) / &
                ( (r_xs(1)-r_xf(1))**2 + (r_xs(2)-r_xf(2))**2 )
           
           !-------------------------------------------------
           ! Ratio of the distances, used for extrapolation.
           !-------------------------------------------------
           
           corr = d_u - 1.0_MK
           
        END IF
        
        !----------------------------------------------------
        ! Set the maximum ratio for the extrapolation.
        !----------------------------------------------------
        
        IF ( corr > mcf_colloid_dist_ratio ) THEN
           
           corr = mcf_colloid_dist_ratio
           
        ELSE IF ( corr < -mcf_colloid_dist_ratio ) THEN
           
           corr = -mcf_colloid_dist_ratio
           
        END IF
        
        !----------------------------------------------------
        ! Get equavilent translational velocity of the
        ! surface point s, if it is rotating.
        !----------------------------------------------------
        
        IF ( this%rotate ) THEN
           
           CALL tool_cross_product(this%tool,&
                this%omega(1:3,sid_c), r_xs(1:3),&
                r_vs(1:3),stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Extrapolated velocity for the colloid boundary 
        ! particle, considering the movment of colloid.
        !----------------------------------------------------
        
        vc(1:2) = -corr *&
             (vf(1:2)-vcoll(1:2)-r_vs(1:2)) + &
             vcoll(1:2) + r_vs(1:2)
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_noslip_Morris_star
      
