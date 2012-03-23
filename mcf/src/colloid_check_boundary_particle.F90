      SUBROUTINE colloid_check_boundary_particle(this,&
           p_x,l_sur,l_out,l_in,c_sid,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_check_boundary_particle
        !----------------------------------------------------
        !
        ! Purpose     : check if a particle is inside a
        !               colloid geometry or 
        !               surrounding a colloid surface;
        !               
        !
        ! Reference   : ellipse/ellipsoid in Wikipedia.
        !
        ! Remark      : 1 :
        !               In case of periodic or Lees-Edwards
        !               boundaries, the images of colloid's
        !               center have to be taken into account
        !               to decide if a potential SDPD boundary
        !               particle is inside the geometry of
        !               a colloid.
        !               For 2D, maximum 3**2=9->8 images;
        !               For 3D, maximum 3**3=27->26 images.
        !              
        !               2 :
        !               In order to prevent the fluid
        !               particle being too close to surface
        !               of a colloid initialy,
        !               'dout' restricts distance
        !               for the fluid particle from surface.
        !
        !               3: 
        !               if a boundary particle is far from
        !               the surface(din or cut_off), 
        !               it is computationally useless.
        !               However, the mass distribution is
        !               still uniform, as if there are
        !               boundary particles far from the surface.
        !
        !
        ! Revisions   : V0.3 14.3.2012, make void for 3D ellipsoid
        !               (prolate and oblate) for computatonal
        !               efficiency.
        !
        !               V0.2 8.3.2012, refine the dicolloid part.
        !
        !               V0.1 10.12 2009, original version.
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
        ! this      : object of colloid class.
        ! p_x       : potential boundary particle's position.
        !
        ! Output  
        !
        ! l_sur     : particle is inside geometry surface.
        ! l_out     : particle is inside geometry surface+dout.
        ! l_in      : particle is inside inner ring, which
        !             is(probaly more than) cut off far 
        !             from surface.
        ! c_sid     : colloid's ID.
        !
        ! stat_info : status of this routine.
        !----------------------------------------------------
        
        TYPE(colloid), INTENT(IN)               :: this
        REAL(MK),DIMENSION(:),INTENT(IN)        :: p_x
        LOGICAL, INTENT(OUT)                    :: l_sur
        LOGICAL, INTENT(OUT)                    :: l_out
        LOGICAL, INTENT(OUT)                    :: l_in
        INTEGER, INTENT(OUT)                    :: c_sid
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here.
        !
        ! rp_x : relative position vector to center of
        !        the colloid.
        ! d_pc : distance of potential boundary particle to
        !        the center of the colloid.
        ! d_sc : distance of from surface to the center of
        !        the colloid.
        ! d_ps : shortest distance of particle to the surface.
        ! theta: the angel between vertor from particle to
        !        center of colloid and x+ direction.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        REAL(MK), DIMENSION(3)                  :: x_coll
        REAL(MK), DIMENSION(3)                  :: rp_x
        REAL(MK), DIMENSION(3)                  :: v_coll
        
        REAL(MK)                                :: a,b,c,d
        REAL(MK)                                :: d_pc
        REAL(MK)                                :: d_pc1,d_pc2
        REAL(MK)                                :: d_sc
        REAL(MK)                                :: d_ps
        REAL(MK)                                :: theta, phi
        REAL(MK)                                :: rp_xy
        REAL(MK), DIMENSION(3)                  :: s_x
        REAL(MK), DIMENSION(3)                  :: center1, center2
        INTEGER                                 :: i,j
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     =  0
        stat_info_sub =  0
        
        l_sur    = .FALSE.
        l_out    = .FALSE.
        l_in     = .FALSE.
        c_sid    = 0
        dim      = this%num_dim
        rp_x(:)  = 0.0_MK

        
        DO i = 1, this%num_colloid
           
           !-------------------------------------------------
           ! Consider the relative positions of nearest 
           ! image centers, including box(cell) itself.
           ! Get relative displacement of potential boundary
           ! particle to the center of nearest image.
           !-------------------------------------------------
           
           CALL colloid_nearest_image(this,p_x(1:dim),i,&
                x_coll(1:dim),rp_x(1:dim),v_coll(1:dim),&
                stat_info_sub)
           
           !-------------------------------------------------
           ! According to different shapes, determine particle
           ! p's location.
           !-------------------------------------------------
           
           a = this%radius(1,i)
           b = this%radius(2,i)
           
           SELECT CASE ( this%shape(i) )
              
           CASE ( mcf_colloid_shape_cylinder )
                 
              !----------------------------------------------
              ! cylinder 2D/3D
              !----------------------------------------------
              
              !----------------------------------------------
              ! For now, assuming it is along z direction
              ! periodically.
              ! Calculate the distance of p from the center
              ! of a colloid, or its image.
              ! Then calculate the angle from x+ direction.
              !----------------------------------------------
           
              d_pc = SQRT(DOT_PRODUCT(rp_x(1:2),rp_x(1:2)))
              
              IF ( d_pc <= a + this%dout ) THEN
                 
                 l_out = .TRUE.
                 
                 IF ( d_pc <= a ) THEN
                    
                    l_sur = .TRUE.
                    
                    IF ( d_pc <= a - this%din ) THEN
                       
                       l_in = .TRUE.
                       
                    END IF ! d_pc <= ra - din
                    
                 END IF ! d_pc <= ra
                 
              END IF ! d_pc <= ra - dout

           CASE ( mcf_colloid_shape_sphere )
                 
              !----------------------------------------------
              ! disk 2D/sphere 3D
              !----------------------------------------------
              !----------------------------------------------
              ! Calculate the distance of p from the center
              ! of a colloid, or its image.
              ! Then calculate the angle from x+ direction.
              !----------------------------------------------
           
              d_pc = SQRT(DOT_PRODUCT(rp_x(1:dim),rp_x(1:dim)))
           
              
              IF ( d_pc <= a + this%dout ) THEN
                 
                 l_out = .TRUE.
                 
                 IF ( d_pc <= a ) THEN
                    
                    l_sur = .TRUE.
                    
                    IF ( d_pc <= a - this%din ) THEN
                       
                       l_in = .TRUE.
                       
                    END IF ! d_pc <= ra - din
                    
                 END IF ! d_pc <= ra
                 
              END IF ! d_pc <= ra - dout
              
              
           CASE ( mcf_colloid_shape_ellipsoid )
              
              !----------------------------------------------
              ! ellipse 2D/ellipsoid 3D
              !----------------------------------------------
              !----------------------------------------------
              ! Calculate the distance of p from the center
              ! of a colloid, or its image.
              ! Then calculate the angle from x+ direction.
              !----------------------------------------------
           
              d_pc = SQRT(DOT_PRODUCT(rp_x(1:dim),rp_x(1:dim)))
           
              IF ( dim == 2 ) THEN
                    
                 !-------------------------------------------
                 ! 2D ellipse:
                 ! get angle in polar coordinate of 
                 ! the potential particle.
                 !-------------------------------------------
                 
                 theta = colloid_polar_angle(rp_x(1),rp_x(2))
                         
                 !----------------------------------------------
                 ! Get distance of the point at angle theta
                 ! on the surface to the center.
                 !----------------------------------------------
                 
                 d_sc = colloid_polar_ellipse_r(a,b, &
                      theta,this%acc_vector(4,i))
                 
                 !----------------------------------------------
                 ! a particle closer than the surface point
                 ! is inside ellipse.
                 !----------------------------------------------
                 
                 IF ( d_pc <= d_sc ) THEN
                    
                    l_sur = .TRUE.
                    
                    !----------------------------------------
                    ! Get distance of the point at angle theta
                    ! on the inner ring surface to the center,
                    ! and the point is saved in s_x(:).
                    !----------------------------------------
                 
                    CALL colloid_cartesian_ellipse_shortestD(&
                         a,b,this%acc_vector(4,i), &
                         rp_x(1),rp_x(2),s_x(1),s_x(2),&
                         d_ps,stat_info_sub)
                    
                    IF ( stat_info_sub /=0 ) THEN
                       PRINT *, __FILE__, ":", __LINE__
                       stat_info = -1
                       GOTO 9999
                    END IF
                    
                    !-------------------------------------------
                    ! particle closer than the inner ring
                    ! (din distance away from surface inside)
                    ! surface point is considered 
                    ! compuationally useless.
                    !-------------------------------------------
                    
                    IF ( d_ps > this%din ) THEN
                       
                       l_in = .TRUE.
                       
                    END IF
                    
                 END IF  ! l_sur
                 
              ELSE IF ( dim == 3 ) THEN
                 
                 c = this%radius(3,i)
                 
                 !-------------------------------------------
                 ! Transpose it to the first orientation.
                 !-------------------------------------------
                 
                 rp_x(1:dim) = MATMUL(&
                      TRANSPOSE(this%acc_matrix(1:dim,1:dim,i)),&
                      rp_x(1:dim) )
                 
                 !-------------------------------------------
                 ! Determine the particles's location.
                 !-------------------------------------------
                 
                 rp_xy = SQRT(rp_x(1)**2+rp_x(2)**2)
                 theta = colloid_polar_angle(rp_x(1),rp_x(2))
                 phi   = colloid_polar_angle(rp_x(3),rp_xy)
                 d_sc  = &
                      colloid_spherical_ellipsoid_r(a,b,c, &
                      theta,phi)
                 
                 IF ( d_pc <= d_sc ) THEN
                    
                    l_sur = .TRUE.
                    
                    !----------------------------------------
                    ! Get distance of the point at angle theta
                    ! on the inner ring surface to the center,
                    ! and the point is saved in s_x(:).
                    !----------------------------------------
                    
                    CALL colloid_cartesian_ellipsoid_shortestD(&
                         this,a,b,c,rp_x(1:3),s_x(1:3),&
                         d_ps,stat_info_sub)
                    
                    IF ( stat_info_sub /=0 ) THEN
                       PRINT *, __FILE__, ":", __LINE__
                       stat_info = -1
                       GOTO 9999
                    END IF
                    
                    !-------------------------------------------
                    ! particle closer than the inner ring
                    ! (din distance away from surface inside)
                    ! surface point is considered 
                    ! compuationally useless.
                    !-------------------------------------------
                    
                    IF ( d_ps > this%din ) THEN
                       
                       l_in = .TRUE.
                       
                    END IF
                    
                 END IF ! d_pc <= d_sc
                 
              END IF ! dim
              
           CASE ( mcf_colloid_shape_dicolloid )
              
              !----------------------------------------------
              ! 2D dicolloid with two disks 
              ! 3D dicolloid with two spheres.
              !----------------------------------------------

              d = a - b
              
              !----------------------------------------------
              ! Calculate the distance of p from the center
              ! of a colloid, or its image.
              ! Then calculate the angle from x+ direction.
              !----------------------------------------------
           
              d_pc = SQRT(DOT_PRODUCT(rp_x(1:dim),rp_x(1:dim)))
              
              !----------------------------------------------
              ! Transpose it to the first orientation. 
              !----------------------------------------------
              
              rp_x(1:dim) = MATMUL(&
                   TRANSPOSE(this%acc_matrix(1:dim,1:dim,i)),rp_x(1:dim))
             
              center1(2:dim) = 0.0_MK
              center2(2:dim) = 0.0_MK
              center1(1) = -d
              center2(1) =  d
              
              d_pc1 = &
                   SQRT(DOT_PRODUCT(rp_x(1:dim)-center1(1:dim),&
                   rp_x(1:dim)-center1(1:dim)))
              d_pc2 = &
                   SQRT(DOT_PRODUCT(rp_x(1:dim)-center2(1:dim),&
                   rp_x(1:dim)-center2(1:dim)))
              
              IF ( d_pc1 <= b + this%dout .OR. &
                   d_pc2 <= b + this%dout  ) THEN
                 
                 l_out = .TRUE.
                 
                 IF ( d_pc1 <= b .OR. &
                      d_pc2 <= b ) THEN
                    
                    l_sur = .TRUE.
                    
                    !----------------------------------------
                    ! Check particle further inside of two
                    ! surfaces is computationally useless.
                    !----------------------------------------
                       
                    IF ( d_pc1 <= b - this%din .OR. &
                         d_pc2 <= b - this%din ) THEN
                       
                       l_in = .TRUE.
                       
                    END IF ! d_pc <= ra - din
                    
                    !----------------------------------------
                    ! Further check if the particle is in the
                    ! overlapping region of two spheres and
                    ! far away from joining point(2D) or
                    ! circle(3D),
                    ! it is also computationally useless.
                    !----------------------------------------
                    
                    IF ( .NOT. l_in .AND. &
                         d_pc1 <= b .AND. &
                         d_pc2 <= b ) THEN
                       
                       !-------------------------------------
                       ! first get p's prejection point
                       ! on y axis(2D)/y-z plane(3D), p1.
                       ! then calculate p1 distance from
                       ! the joining point(2D) or circle(3D).
                       ! Then calculate p's distance from
                       ! joining point or circle.
                       !-------------------------------------
                       
                       d_pc = 0.0_MK
                       
                       DO j = 2, dim
                          
                          d_pc = d_pc + rp_x(j)**2
                          
                       END DO

                       d_pc = SQRT(d_pc)
                       
                       d_pc = SQRT(b**2-d**2) - d_pc
                       
                       d_pc = SQRT(rp_x(1)**2+d_pc**2)
                       
                       IF ( d_pc > this%din ) THEN
                          
                          l_in = .TRUE.
                          
                       END IF
                       
                    END IF
                    
                 END IF ! d_pc <= ra
                 
              END IF ! d_pc <= ra - dout
              
           CASE ( mcf_colloid_shape_star )
              
              !----------------------------------------------
              ! Calculate the distance of p from the center
              ! of a colloid, or its image.
              ! Then calculate the angle from x+ direction.
              !----------------------------------------------
           
              d_pc = SQRT(DOT_PRODUCT(rp_x(1:dim),rp_x(1:dim)))
           
              !----------------------------------------------
              ! 2D star with different frequency.
              !
              ! At angel theta, caculate the point on the
              ! surface to the center.
              !
              ! 3D assuming rotating 2D shape wiht x-axis for
              ! 2pi.
              !----------------------------------------------
              
              theta = colloid_polar_angle(rp_x(1),rp_x(2))
              
              d_sc  = colloid_polar_star_r(&
                   a,b,&
                   REAL(this%freq(i),MK),theta,this%acc_vector(4,i))
              
              !----------------------------------------------
              ! If the particle is closer than the point on 
              ! the surface at same angel theta, it is inside.
              !----------------------------------------------
              
              IF ( d_pc <= d_sc ) THEN
                 
                 l_sur = .TRUE.
                 
                 !-------------------------------------------
                 ! At angel theta, caculate the point on the
                 ! inner ring surface to the center.
                 !-------------------------------------------
                 
                 CALL colloid_polar_star_shortestD(&
                      this%radius(1,i),this%radius(2,i),&
                      REAL(this%freq(i),MK),&
                      this%acc_vector(4,i),rp_x(1),rp_x(2),&
                      s_x(1),s_x(2),d_ps,stat_info_sub)
               
                 IF ( stat_info_sub /= 0 ) THEN
                    PRINT *, "colloid_check_boundary_particle : ", &
                         "Finding star shortest D wrong !"
                    stat_info = -1
                    GOTO 9999                    
                 END IF
                 
                 !-------------------------------------------
                 ! particle closer than the inner ring
                 ! (din far away from surface)
                 ! surface point is considered 
                 ! compuationally useless.
                 !-------------------------------------------
                 
                 IF ( d_ps >= this%din ) THEN
                    
                    l_in = .TRUE.
                    
                 END IF
                 
              END IF  ! l_sur

           CASE DEFAULT
              
              PRINT *, "colloid_check_boundary_particle: ", &
                   "No shape available !"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! shape
           
           !-------------------------------------------------
           ! If the boundary particle is found to be inside
           ! of any surface, stop searching further.
           !-------------------------------------------------
           
           IF ( l_sur ) THEN
              
              c_sid  = i           
              
              EXIT
              
           END IF
           
        END DO ! i = 1, num_colloid
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  colloid_check_boundary_particle
