      SUBROUTINE colloid_noslip_Morris_ellipsoid(this, &
           xf,xc,vf,vc,sid_c,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_noslip_Morris_ellipsoid
        !----------------------------------------------------
        ! Purpose     : For 3D ellipsoid.
        !
        ! Remark      : As the ellipsoid can rotate,
        !               quantities, such as, position,
        !               translational/rotational velocity
        !               have to be transfered to body-attached
        !               coordinate A.
        !               After the extrapolated velocity
        !               is calculated in A, it has to be
        !               transfered back to fixed coordinate.
        !
        ! Revision    : V0.1 14.03.2012, orignal version.
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
        ! dim     : number of dimension.
        ! xcoll   : position of the center of the colloid.
        ! vcoll   : velocity of the center of the colloid.
        ! r_xf    : relative position from fluid particle 
        !           to the center.
        ! r_xc    : relative position from colloid boundary
        !           particle to center.
        ! r_vc    : relative velocity of colloid boundary 
        !           particle to the center.
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
        REAL(MK), DIMENSION(3)                  :: ocoll
        REAL(MK), DIMENSION(3)                  :: r_xf
        REAL(MK), DIMENSION(3)                  :: r_xc
        REAL(MK), DIMENSION(3)                  :: r_vc
        REAL(MK), DIMENSION(3)                  :: r_xs
        REAL(MK), DIMENSION(3)                  :: r_vs        
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
        !PRINT *, "xc, sid_c, xcoll, r_xc, vcoll: ", &
        !     xc(:), sid_c, xcoll(:), r_xc(:), vcoll(:)
        !----------------------------------------------------
        ! Get relative position of fluid particle to 
        ! the center of the geometry.
        !----------------------------------------------------
        
        r_xf(1:dim) = xf(1:dim) - xcoll(1:dim)
        !PRINT *, "r_xf: ", r_xf(:)
        !----------------------------------------------------
        ! Using body attached coordinate,
        ! the ellipsoid is in its frist orientation,
        ! i.e.,
        ! a along x direction;
        ! b along y direction;
        ! c along z direction.
        !
        ! Get positions of fluid and boundary particle 
        ! in this body attached coordinate.
        !----------------------------------------------------
        
        r_xf(1:dim) = MATMUL(&
             TRANSPOSE(this%acc_matrix(1:dim,1:dim,sid_c)),&
             r_xf(1:dim) )
        r_xc(1:dim) = MATMUL(&
             TRANSPOSE(this%acc_matrix(1:dim,1:dim,sid_c)),&
             r_xc(1:dim) )
        !PRINT *, "acc_matrix: ", this%acc_matrix(:,:,sid_c)
        !PRINT *, "r_xf, r_xc: ", r_xf(:), r_xc(:)
                
        !----------------------------------------------------
        ! Find the shortest distance from the fluid particle
        ! to the surface of the ellipsoid.
        ! r_xs(:) contains the surface point.
        !----------------------------------------------------
        
        CALL colloid_cartesian_ellipsoid_shortestD(&
             this, &
             this%radius(1,sid_c),&
             this%radius(2,sid_c),&
             this%radius(3,sid_c),&
             r_xf(1:3),r_xs(1:3),&
             d_fs,stat_info_sub)
        
        IF ( stat_info_sub /=0 ) THEN
           PRINT *, __FILE__, ":", __LINE__
           stat_info = -1
           GOTO 9999
        END IF
        !PRINT *, "r_xf, r_xs, d_fs: ", r_xf(:), r_xs(:), d_fs
        
#if 0
        !----------------------------------------------------
        ! Check if the fluid particle is inside of 
        ! the geometry.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Check if the boundary particle is outside of
        ! the geometry.
        !----------------------------------------------------
        
#endif
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
           !PRINT *, "nvector: ", nvector(:)
           !-------------------------------------------------
           ! Map (r_xs-r_xc) onto nvector and 
           ! calculate its length.
           !-------------------------------------------------
        
           d_cs = &
                DOT_PRODUCT((r_xs(1:dim)-r_xc(1:dim)),&
                nvector(1:dim))
           !PRINT *, "d_cs: ", d_cs
           !-------------------------------------------------
           ! Ratio of the distances, used for extrapolation.
           !-------------------------------------------------
           
           corr = d_cs / d_fs
           !PRINT *, "corr: ", corr
           !-------------------------------------------------
           ! Set the maximum ratio of the extrapolation.
           !-------------------------------------------------
           
           IF ( corr > mcf_colloid_dist_ratio_max ) THEN
              
              corr = mcf_colloid_dist_ratio_max
              
           END IF
           
        END IF ! d_fs == 0
        
        !----------------------------------------------------
        ! Calculate the translational/rotational velcoity
        ! of the colloid in body-attached coordinate.
        !----------------------------------------------------
        
        vcoll(1:dim) = MATMUL(&
             TRANSPOSE(this%acc_matrix(1:dim,1:dim,sid_c)),&
             vcoll(1:dim) )
        ocoll(1:dim) = MATMUL(&
             TRANSPOSE(this%acc_matrix(1:dim,1:dim,sid_c)),&
             this%omega(1:dim,sid_c,1) )
        !PRINT *, "vcoll, ocoll: ", vcoll(:), ocoll(:)
        !----------------------------------------------------
        ! Extrapolate velocity for the colloid boundary 
        ! particle, considering the movment of colloid.
        !----------------------------------------------------
        
        r_vs(:) = 0.0_MK
        
        CALL tool_cross_product(this%tool,&
             ocoll(1:dim), r_xs(1:dim),&
             r_vs(1:dim),stat_info_sub)
        !PRINT *, "r_vs: ", r_vs(:)
        vc(1:dim) = -corr *&
             (vf(1:dim)-vcoll(1:dim)-r_vs(1:dim)) + &
             vcoll(1:dim) + r_vs(1:dim)
        !PRINT *, "vc: ", vc(:)
        !----------------------------------------------------
        ! Transfer the result back to the fixed coordinate.
        !----------------------------------------------------
        
        vc(1:dim) = &
             MATMUL(this%acc_matrix(1:dim,1:dim,sid_c),&
             vc(1:dim) )
        !PRINT *, "vc: ", vc(:)
        
9999    CONTINUE
        !STOP
        RETURN
        
      END SUBROUTINE colloid_noslip_Morris_ellipsoid
      
