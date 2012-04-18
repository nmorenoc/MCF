      SUBROUTINE colloid_noslip_Morris_cylinder_2D(this, &
           xf,xc,vf,vc,sid_c,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_noslip_Morris_cylinder
        !----------------------------------------------------
        ! Purpose     : For 2D cylinder velocity
        !               no-slip boundary condition.
        !
        ! Remark      : it is identical to 2D disk.
        !
        ! Revision    : V0.3 4.3.2010, small bug is fixed
        !               according to Bian et al. 2012
        !               paper.
        !               The bug was irrelevant for 
        !               irrotational colloid.
        !               
        !               V0.2 27.11.2009, including
        !               Lees-Edwards boundary.
        !
        !               V0.1 01.03.2009, orignal version,
        !               with periodic boundary.
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
        ! d_ccoll : distance from c to the image center.
        ! nvector : unit vector pointing from the image
        !           center to particle f.
        ! d_cn    : distance mapping d_ccoll to nvector.
        ! d_fs    : distance from f to tangent surface.
        ! d_cs    : distance from c to tangent surface.
        ! corr    : d_cs/d_fs.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        REAL(MK), DIMENSION(3)                  :: xcoll
        REAL(MK), DIMENSION(3)                  :: vcoll
        REAL(MK), DIMENSION(3)                  :: r_xf
        REAL(MK), DIMENSION(3)                  :: r_xc
        REAL(MK), DIMENSION(3)                  :: r_xs
        REAL(MK), DIMENSION(3)                  :: r_vc
        REAL(MK), DIMENSION(3)                  :: r_vs
        REAL(MK)                                :: d_fcoll
        REAL(MK)                                :: d_ccoll
        REAL(MK), DIMENSION(3)                  :: nvector
        REAL(MK)                                :: d_cn
        REAL(MK)                                :: d_fs
        REAL(MK)                                :: d_cs
        REAL(MK)                                :: corr
          
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
        nvector(:) = 0.0_MK

        !----------------------------------------------------
        ! Get the nearest image colloid center's to the
        ! boundary particle.
        !----------------------------------------------------
        
        CALL colloid_in_nearest_image(this,xc(1:dim),sid_c,&
             xcoll(1:dim),r_xc(1:dim),vcoll(1:dim),&
             stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, __FILE__, ":", __LINE__
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Get relative position of fluid particle to the 
        ! image center.
        !----------------------------------------------------
        
        r_xf(1:2) = xf(1:2) - xcoll(1:2)
        
        !----------------------------------------------------
        ! Distance between f and the image center.
        !----------------------------------------------------
        
        d_fcoll = SQRT(DOT_PRODUCT(r_xf(1:dim),r_xf(1:dim)))
        
        !----------------------------------------------------
        ! Check if penetration happens, i.e., a fluid 
        ! particle goes inside a colloid object. If yes, 
        ! simply assign the same velocity as the colloid.
        !----------------------------------------------------
        
        IF( d_fcoll <= this%radius(1,sid_c) ) THEN
           
           vc(1:dim) = vcoll(1:dim)
           
           !----------------------------------------------------
           ! Get equavilent translational velocity of the
           ! angular velocity, if it is rotating.
           !----------------------------------------------------
           
           IF ( this%rotate ) THEN
              
              r_vc(:) = 0.0_MK
              
              CALL tool_cross_product(this%tool,&
                   this%omega(1:3,sid_c,1), r_xc(1:3),&
                   r_vc(1:3),stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, __FILE__, ":", __LINE__
                 stat_info = -1
                 GOTO 9999
              END IF
              
              vc(1:dim) = vc(1:dim) + r_vc(1:dim)
              
           END IF
           
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
        
        IF( d_ccoll > this%radius(1,sid_c) + this%dout ) THEN
           PRINT *, "colloid_noslip_Morris_sphere : ", &
                "Colloid particle goes out of colloid disk/sphere !"
           PRINT *, "dist,radius,xc,xcolloid,xcoll :", &
                d_ccoll,this%radius(1,sid_c),xc(1:dim),&
                this%x(1:dim,sid_c),xcoll(1:dim)
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Normalize r_xf and get nvector.
        !----------------------------------------------------
        
        nvector(1:2) = r_xf(1:2) / d_fcoll
        
        !----------------------------------------------------
        ! Map r_xc on nvector and caculate the length.
        !----------------------------------------------------
        
        d_cn = DOT_PRODUCT(r_xc(1:2),nvector(1:2))
        
        !----------------------------------------------------
        ! Distance from the colloid boundary particle to 
        ! the tangent surface of the colloid.
        !----------------------------------------------------
        
        d_cs = this%radius(1,sid_c) - d_cn
        
        !----------------------------------------------------
        ! Distance from the fluid particle to the 
        ! tangent surface of the colloid.
        ! Calculate also r_xs vector.
        !----------------------------------------------------
        
        d_fs = d_fcoll - this%radius(1,sid_c)
        
        r_xs(1:dim) =this%radius(1,sid_c)/d_fcoll*r_xf(1:dim)
        
        !----------------------------------------------------
        ! If the fuild particle lies exactly on or in
        ! the surface, assign a minimum distance.
        !
        ! f inside surface should not ever happen(d_fs<0).
        !----------------------------------------------------
        
        IF( ABS(d_fs) <  mcf_machine_zero ) THEN
           
           d_fs = mcf_machine_zero
           
           corr = mcf_colloid_dist_ratio_max
           
        ELSE
           
           !-------------------------------------------------
           ! Ratio of the distances, used for extrapolation.
           !-------------------------------------------------
           
           corr = d_cs / d_fs
           
           !----------------------------------------------------
           ! Set the maximum ratio for the extrapolation.
           !----------------------------------------------------
           
           IF ( corr > mcf_colloid_dist_ratio_max ) THEN
              
              corr = mcf_colloid_dist_ratio_max
              
           ELSE IF ( corr < -mcf_colloid_dist_ratio_max ) THEN
              
              corr = -mcf_colloid_dist_ratio_max
              
           END IF

        END IF
        
        !----------------------------------------------------
        ! Get equavilent translational velocity of the
        ! surface point s, if it is rotating.
        !----------------------------------------------------
        
        IF ( this%rotate ) THEN
           
           CALL tool_cross_product(this%tool,&
                this%omega(1:3,sid_c,1), r_xs(1:3),&
                r_vs(1:3),stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Extrapolated velocity for the colloid boundary 
        ! particle, considering the movment of colloid.
        !----------------------------------------------------
        
        vc(1:dim) = -corr *&
             (vf(1:dim)-vcoll(1:dim)-r_vs(1:dim)) + &
             vcoll(1:dim) + r_vs(1:dim)
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_noslip_Morris_cylinder_2D
