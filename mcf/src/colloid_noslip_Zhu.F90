      SUBROUTINE colloid_noslip_Zhu(this,xf,xc,vf,vc,sid_c,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_noslip_Zhum
        !----------------------------------------------------
        ! Purpose     : Private subroutine implementing
        !               the noslip condition from 
        !               Zhu et al 1999.
        !               A Spheric 3D/ cylinder 2D 
        !               colloid object.
        !
        ! Revision    : V0.2 27.11.2009, including
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
        ! Arguments
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
        INTEGER                                 :: dim,i,k
        INTEGER                                 :: num_le
        REAL(MK), DIMENSION(3)                  :: x_coll
        REAL(MK), DIMENSION(3)                  :: v_coll
        REAL(MK),DIMENSION(:,:),POINTER         :: shear_length
        REAL(MK),DIMENSION(:,:),POINTER         :: shear_v
        REAL(MK), DIMENSION(3)                  :: length
        REAL(MK)                                :: fcoll_dist
        REAL(MK)                                :: ccoll_dist
        REAL(MK), DIMENSION(3)                  :: nvector
        REAL(MK)                                :: cn_dist
        REAL(MK)                                :: fsur_dist,csur_dist
        REAL(MK)                                :: corr        
        REAL(MK), DIMENSION(3)                  :: r
        REAL(MK), DIMENSION(3)                  :: v_r

        !------------------------------------------
        ! Initialization.
        !------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0 

        dim = this%num_dim
        
        NULLIFY(shear_length)
        NULLIFY(shear_v)
        
        !------------------------------------------
        ! If there is Lees-Edwards boundary,
        ! get its shear length and velocity.
        !------------------------------------------

        num_le = boundary_get_num_le( this%boundary,&
             stat_info_sub)
        
        IF ( num_le > 0 ) THEN
           
           CALL boundary_get_shear_length(this%boundary,&
                shear_length,stat_info_sub)
           CALL boundary_get_shear_v(this%boundary,&
                shear_v,stat_info_sub)

        END IF
        
        !------------------------------------------
        ! Get the length of computational box.
        !------------------------------------------
        
        length(1:dim) = &
             this%max_phys(1:dim) - this%min_phys(1:dim)
        
        !------------------------------------------
        ! Get original collloid center's position
        ! and velocity.
        !------------------------------------------
        
        x_coll(1:dim) = this%x(1:dim,sid_c)
        v_coll(1:dim) = this%v(1:dim,sid_c,1)
        
        !---------------------------------------
        ! In context of periodic or Lees-Edwards
        ! boundaries, we have to consider the 
        ! colloid boundary particle which is not 
        ! only ghost particle but also periodic 
        ! images in order to get correct positions.
        !---------------------------------------
        
        Do i = 1, dim
           
           !---------------------------------------
           ! Check if the colloid boundary particle 
           ! is in around the colloid center.
           ! If yes, use the center coordinate of 
           ! the colloid.
           ! If no, translate the center coodinate
           ! of the colloid to get the image center.
           !---------------------------------------
           
           IF( xc(i) - x_coll(i) > &
                this%radius(1,sid_c) + this%dout) THEN
              
              SELECT CASE ( this%bcdef(2*i-1) ) 
                 
              CASE ( ppm_param_bcdef_periodic )
                 
                 x_coll(i) = x_coll(i) + &
                      length(i)
                 
              CASE ( ppm_param_bcdef_LE )
                 
                 DO k = 1, dim
                    
                    IF ( k == i ) THEN
                       
                       x_coll(k) = x_coll(k) + length(k)
                       
                    ELSE
                       
                       x_coll(k) = &
                            MODULO(x_coll(k) - shear_length(k,2*i-1), &
                            length(k) )
                       
                       v_coll(k) = v_coll(k) + &
                            ( shear_v(k,2*i) - shear_v(k,2*i-1) )
                       
                    END IF
                    
                 END DO ! k
                 
              END SELECT ! bcdef
              
           ELSE IF( x_coll(i) - xc(i) >&
                this%radius(1,sid_c) + this%dout ) THEN
              
              SELECT CASE ( this%bcdef(2*i) )
                 
              CASE ( ppm_param_bcdef_periodic )
                 
                 x_coll(i) = x_coll(i) - &
                      length(i)
                 
              CASE ( ppm_param_bcdef_LE )
                 
                 DO k = 1, dim
                    
                    IF ( k == i ) THEN
                       
                       x_coll(k) = x_coll(k) - length(k)
                       
                    ELSE
                       
                       x_coll(k) = &
                            MODULO(x_coll(k) - shear_length(k,2*i), &
                            length(k) )
                       
                       v_coll(k) = v_coll(k) + &
                            ( shear_v(k,2*i-1) - shear_v(k,2*i) )
                       
                    END IF
                    
                 END DO
                 
              END SELECT
              
           END IF
           
        END DO
        
        !---------------------------------------
        ! Distance between a fluid particle and 
        ! the center of a colloid.
        !--------------------------------------- 
        
        fcoll_dist = 0.0_MK
        DO i = 1, dim
           fcoll_dist = fcoll_dist + ( xf(i) - x_coll(i) )**2
        END DO
        fcoll_dist = SQRT(fcoll_dist)
        
        !-----------------------------------------------
        ! Check if penetration happens, i.e.,
        ! a fluid particle goes inside a colloid object.
        ! If yes, simply assign the same velocity as
        ! the colloid.
        !-----------------------------------------------
        
        IF(fcoll_dist <=  this%radius(1,sid_c) ) THEN
           vc(1:dim) = v_coll(1:dim)
           !PRINT *, "Fluid particle inside the colloid sphere!"          
           !stat_info = -1
           !GOTO 9999           
        END IF
        
        !-----------------------------------------------
        ! Distance between colloid boundary particle 
        ! and the center of a colloid.
        !-----------------------------------------------
        
        ccoll_dist = 0.0_MK    
        DO i =1, dim
           ccoll_dist =  ccoll_dist+ ( xc(i)-x_coll(i) )**2
        END DO
        ccoll_dist = SQRT(ccoll_dist)

        !----------------------------------------------------
        ! Inconsistent movement happens for a boundary 
        ! particle inside a colloid, it is wrong.
        !----------------------------------------------------
        
        IF( ccoll_dist > this%radius(1,sid_c) + this%dout) THEN
           PRINT *, "colloid_noslip_Zhu_sphere : ", &
                "Colloid particle goes out of colloid sphere !"
           PRINT *, "dist,radius, xc,xcoll, x_coll :", &
                ccoll_dist,this%radius(1,sid_c), &           
                xc(1:dim),&
                this%x(1:dim,sid_c),x_coll(1:dim)
           stat_info = -1
           GOTO 9999      
        END IF
        
        !---------------------------------------
        !  Normal vector pointing from center of 
        !  the colloid to the fluid particle.
        !---------------------------------------
        
        nvector(1:dim) = &
             ( xf(1:dim)-x_coll(1:dim) ) /  fcoll_dist  
        
        !------------------------------------------
        ! Vector pointing from centre of the 
        ! colloid to the colloid boundary particle,
        ! then mapped it onto normal vector.
        !------------------------------------------

        cn_dist = 0.0_MK
        DO i =1, this%num_dim           
           cn_dist = cn_dist+ (xc(i)-x_coll(i)) * nvector(i)
        END DO

        !---------------------------------------
        ! Distance from the colloid boundary 
        ! particle to the tangent surface of the 
        ! colloid.
        !---------------------------------------
        
        csur_dist = this%radius(1,sid_c) - cn_dist
        
        !------------------------------------------
        ! Distance from the fluid particle to the 
        ! tangent surface of the colloid.
        !------------------------------------------
        
        fsur_dist = fcoll_dist - this%radius(1,sid_c)
        
        !------------------------------------------
        ! If the fuild particle lies exactly on
        ! the surface, assign a minimum distance.
        !------------------------------------------
        
        IF( fsur_dist <= 0.0_MK) THEN
           fsur_dist = mcf_machine_zero
        END IF
        
        
        !temp = SQRT(3.0_MK) * this%dx / 4.0_MK
        
        IF (fsur_dist  <  this%dout) THEN
           fsur_dist = this%dout
        END IF
        

        !--------------------------------
        ! Ratio of the distances, 
        ! used for extrapolation.
        !--------------------------------
        
        corr = csur_dist/fsur_dist
        

        !------------------------------------------
        ! Get the translational angular velocity
        ! of the colloid, if it is rotating.
        !------------------------------------------

        v_r(1:3) = 0.0_MK
        
        IF ( this%rotate ) THEN
        
           r(1:3) = 0.0_MK
           r(1:dim) = xc(1:dim) - x_coll(1:dim)
           
           CALL tool_cross_product(this%tool,&
                this%omega(1:3,sid_c,1), r(1:3),&
                v_r(1:3),stat_info_sub)

        END IF
        
        !------------------------------------------
        ! Extrapolated velocity for the colloid 
        ! boundary particle, considering the movment 
        ! of the colloid itself.
        !------------------------------------------
        
        vc(1:dim) = -corr *&
             (vf(1:dim)-v_coll(1:dim)-v_r(1:dim)) + &
             v_coll(1:dim) + v_r(1:dim)
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(shear_length)) THEN
           DEALLOCATE(shear_length) 
        END IF
        
        IF(ASSOCIATED(shear_v)) THEN
           DEALLOCATE(shear_v) 
        END IF
        
        RETURN
        
      END SUBROUTINE colloid_noslip_Zhu
      
      
