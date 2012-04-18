      SUBROUTINE colloid_in_nearest_image(this,x,sid, &
        x_image,rx,v_image,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_in_nearest_image
        !----------------------------------------------------
        ! Purpose     : Compute a boundary particle's(inside
        !               the colloid) the nearest image 
        !               of the colloid center.
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     : According to different boundary
        !               conditions, images of colloid have
        !               to be considered.
        !
        ! Revisions   : V0.1 11.10.2010, original version
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
        ! x     : position of a fluid particle.
        ! sid   : species ID of a colloid boundary particle.
        !
        ! Output
        !
        ! x_image : position of nearest image of the colloid
        !           center to the boundary particle.
        ! rx      : relative postion.
        ! v_image : velocity of the nearest image.
        !         
        ! stat_info : status of the routine.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: x
        INTEGER, INTENT(IN)                     :: sid
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: x_image
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: rx
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: v_image
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local parameters
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim, i, k
        REAL(MK)                                :: d_max
        REAL(MK), DIMENSION(3)                  :: length
        INTEGER                                 :: num_le
        REAL(MK), DIMENSION(:,:), POINTER       :: shear_length
        REAL(MK), DIMENSION(:,:), POINTER       :: shear_v
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        dim           = this%num_dim
        
        length(1:dim) = &
             this%max_phys(1:dim) - this%min_phys(1:dim)
        
        num_le        = &
             boundary_get_num_le(this%boundary,stat_info_sub)
        
        NULLIFY(shear_length)
        NULLIFY(shear_v)

        CALL boundary_get_shear_length(this%boundary,&
             shear_length,stat_info_sub)
        CALL boundary_get_shear_v(this%boundary,&
             shear_v,stat_info_sub)
       
        
        !----------------------------------------------------
        ! Calculate maximum possible boundary particle 
        ! displacement in all direction relative to the 
        ! colloid center.
        !----------------------------------------------------

        IF ( dim == 2 ) THEN
           
           SELECT CASE ( this%shape(sid) ) 
              
           CASE ( mcf_colloid_shape_cylinder )
              
              d_max = this%radius(1,sid)
              
           CASE ( mcf_colloid_shape_disk ) 
              
              d_max = this%radius(1,sid)
           
           CASE ( mcf_colloid_shape_ellipse )
              
              d_max = this%radius(1,sid)
           
           CASE (mcf_colloid_shape_dicolloid)
           
              d_max = this%radius(1,sid) 

           CASE (mcf_colloid_shape_star) 
              
              d_max = this%radius(1,sid) + this%radius(2,sid)
              
           CASE DEFAULT
           
              PRINT *, __FILE__, ":", __LINE__, &
                   "No such shape available !"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! shape
           
        ELSE IF ( dim == 3 ) THEN
           
           SELECT CASE ( this%shape(sid) ) 
              
           CASE ( mcf_colloid_shape_cylinder )
              
              d_max = this%radius(1,sid)
              
              IF ( this%radius(3,sid) > this%radius(1,sid) ) THEN
                 
                 d_max = this%radius(3,sid)
                 
              END IF
              
           CASE ( mcf_colloid_shape_sphere ) 
              
              d_max = this%radius(1,sid)
              
           CASE ( mcf_colloid_shape_ellipsoid )
              
              d_max = this%radius(1,sid)
              
           CASE (mcf_colloid_shape_dicolloid)
              
              d_max = this%radius(1,sid) 
              
           CASE DEFAULT
              
              PRINT *, __FILE__, ":", __LINE__, &
                   "No such shape available !"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! shape

        END IF ! dim 
        
        !----------------------------------------------------
        ! Give a tolerance dout for inconsistent movement
        ! of a rigid colloid.
        !----------------------------------------------------
        
        d_max = d_max + this%dout
        
        !----------------------------------------------------
        ! Calculate absolute relative position rx(a guess).
        !----------------------------------------------------
        
        x_image(1:dim) = this%x(1:dim,sid)
        v_image(1:dim) = this%v(1:dim,sid,1)
        rx(1:dim)      = x(1:dim) - x_image(1:dim)
        
        !----------------------------------------------------
        ! If there is Lees-Edwards boundary condition,
        ! handle it first.
        ! Check each dimension for LE boundary condition
        ! and then calculate the relative displacement of
        ! boundary particle to the center of colloid 
        ! at inital position, i.e., not sheared yet.
        !----------------------------------------------------
        
        IF ( num_le > 0 ) THEN
           
           !-------------------------------------------------
           ! Check each direction.
           !-------------------------------------------------
           
           DO i = 1, dim
              
              IF ( rx(i) > d_max .AND. &
                   this%bcdef(2*i-1) == ppm_param_bcdef_LE ) THEN
                 
                 !-------------------------------------------
                 ! boundary particle is far away 
                 ! in upper direction, and lower boundary is
                 ! LE boundary condition, that means
                 ! the boundary particle cross lower
                 ! boundary, but the center of colloid not.
                 ! Adjust colloid center with length(i) in LE 
                 ! direction and with shear_length in
                 ! shear directions. Adjust colloid's velcoity
                 ! according to shear velocity.
                 !-------------------------------------------
                 
                 DO k = 1, dim
                    
                    IF ( k == i ) THEN
                       
                       x_image(k) = x_image(k) + length(k)
                       
                    ELSE
                       
                       x_image(k) = &
                            MODULO( x_image(k) + shear_length(k,2*i), &
                            length(k) )
                       
                       v_image(k) = v_image(k) + &
                            ( shear_v(k,2*i) - shear_v(k,2*i-1) )
                       
                    END IF ! k == i
                    
                 END DO ! k = 1, dim
                 
                 EXIT ! Once found, ignore other dimensions.
                 
              ELSE IF ( -rx(i) > d_max .AND. &
                   this%bcdef(2*i) == ppm_param_bcdef_LE ) THEN
                 
                 !----------------------------------------------
                 ! boundary particle is far away 
                 ! in lower direction, and upper boundary is
                 ! LE boundary condition, that means
                 ! the boundary particle cross upper
                 ! boundary, but the center of colloid not.
                 ! Adjust colloid center with -length(i) in LE 
                 ! direction and with shear_length in
                 ! shear directions. Adjust colloid's velocity
                 ! according to shear velocity.
                 !----------------------------------------------
                 
                 DO k = 1, dim
                    
                    IF ( k == i ) THEN
                       
                       x_image(k) = x_image(k) - length(k)
                       
                    ELSE
                       
                       x_image(k) = &
                            MODULO( x_image(k) + shear_length(k,2*i-1), &
                            length(k) )
                       
                       v_image(k) = v_image(k) + &
                            ( shear_v(k,2*i-1) - shear_v(k,2*i) )
                       
                    END IF ! k == i
                    
                 END DO !  k = 1, dim
                 
                 EXIT ! once found, ignore other dimensions.
                 
              END IF ! rx(i) > d_max or -rx(i) > d_max
              
           END DO ! i = 1, num
           
           !-------------------------------------------------
           ! After considering LE boundary condition with 
           ! sheared length, calculate new relative position.
           !-------------------------------------------------
           
           rx(1:dim)    = x(1:dim) - x_image(1:dim)
           
        END IF ! num_le > 0
        
        
        !----------------------------------------------------
        ! Check each dimension for too big displacement of
        ! boundary particle.
        !----------------------------------------------------
        
        DO i = 1, dim
           
           IF ( rx(i) > d_max ) THEN
              
              !----------------------------------------------
              ! boundary particle is far away 
              ! in upper direction, and lower boundary is
              ! periodic boundary condition, that means
              ! the boundary particle cross lower
              ! boundary, but the center of colloid not.
              ! Adjust colloid center with length(i).
              !----------------------------------------------
              
              IF ( this%bcdef(2*i-1) == ppm_param_bcdef_periodic ) THEN
                 
                 x_image(i) = x_image(i) + length(i)
                 
              ELSE
                 
                 PRINT *, "colloid_in_nearest_image : ", &
                      "Can not be non-periodic boundary,", &
                      x(1:dim), sid
                 
                 stat_info = -1
                 GOTO 9999
                 
              END IF ! bcdef(2*i-1)
              
           ELSE IF ( -rx(i) > d_max  ) THEN
              
              !----------------------------------------------
              ! boundary particle is far away 
              ! in lower direction, and upper boundary is
              ! periodic boundary condition, that means
              ! the boundary particle cross upper
              ! boundary, but the center of colloid not.
              ! Adjust colloid center with length(i).
              !----------------------------------------------
             
              IF ( this%bcdef(2*i) ==  ppm_param_bcdef_periodic ) THEN
                 
                 x_image(i) = x_image(i) - length(i)
                 
              ELSE
                 
                 PRINT *,  "colloid_in_nearest_image : ", &
                      "Can not be non-periodic boundary,", &
                      x(1:dim), sid
                 
                 stat_info = -1
                 GOTO 9999
                 
              END IF  ! bcdef(2*i)
              
           END IF ! rx(i) > d_max or -rx(i) > d_max
           
        END DO ! i = 1, dim
        
        !----------------------------------------------------
        ! The actual relative position of the boundary
        ! particle to the nearest image of colloid center
        ! now can be calculated.
        !----------------------------------------------------
        
        rx(1:dim) = x(1:dim) - x_image(1:dim)
        
        
9999    CONTINUE
        
        IF ( ASSOCIATED(shear_length) ) THEN
           DEALLOCATE(shear_length)
        END IF
        
        IF ( ASSOCIATED(shear_v) ) THEN
           DEALLOCATE(shear_v)
        END IF
        
        RETURN
        
      END SUBROUTINE colloid_in_nearest_image
      
      
