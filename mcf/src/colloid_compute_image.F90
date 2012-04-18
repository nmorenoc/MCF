      SUBROUTINE colloid_compute_image(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_image
        !----------------------------------------------------
        !
        ! Purpose     : According to different boundary
        !               conditions, compute the images of
        !               all colloids and save them in the 
        !               object data structure for future usage.
        !               
        !
        ! Reference   : 
        !
        ! Remark      : In case of periodic or Lees-Edwards
        !               boundaries, the images of colloid's
        !               center has to be taken into account.
        !               
        !               For 2D, maximum 3**2=9=8 images;
        !               For 3D, maximum 3**3=27=26 images
        !               (including itself).
        !
        ! Revisions   : V0.1 12.10 2010, original version
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
        !
        ! Output
        !
        ! stat_info : status of this routine.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(OUT)            :: stat_info
                
        !----------------------------------------------------
        ! Local variables
        !
        ! dim    : nubmer of dimension
        ! image  : relative position of image to 
        !          the center box.
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: dim
        INTEGER                         :: num_colloid
        INTEGER                         :: num_le
        INTEGER                         :: num,num_image
        INTEGER                         :: i,j,k
        
        REAL(MK),DIMENSION(3)           :: length
        REAL(MK),DIMENSION(:,:),POINTER :: shear_length
        REAL(MK),DIMENSION(:,:),POINTER :: shear_v
       
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------

        stat_info     = 0
        stat_info_sub = 0
        
        dim         = this%num_dim
        num_colloid = this%num_colloid
        
        !----------------------------------------------------
        ! length : length of the box
        !----------------------------------------------------
        
        length(1:dim) = &
             this%max_phys(1:dim) - this%min_phys(1:dim)
        
        num_le = &
             boundary_get_num_le(this%boundary,stat_info_sub)
        NULLIFY(shear_length)
        NULLIFY(shear_v)
        
        !----------------------------------------------------
        ! Reset image to zero.
        !----------------------------------------------------
        
        this%x_image(:,:,:) = 0.0_MK
        this%v_image(:,:,:) = 0.0_MK
        
        !----------------------------------------------------
        ! If there is L.E. b.c. get shear length and velocity.
        !----------------------------------------------------
        
        IF ( num_le > 0 ) THEN
           
           CALL boundary_get_shear_length(this%boundary,&
                shear_length, stat_info_sub)
           CALL boundary_get_shear_v(this%boundary,&
                shear_v, stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Record the center box first, its relative position
        ! and velocity are zero.
        !----------------------------------------------------
        
        num = 1
        this%x_image(1:dim,1:num_colloid,1) = &
             this%x(1:dim,1:num_colloid)
        this%v_image(1:dim,1:num_colloid,1) = &
             this%v(1:dim,1:num_colloid,1)
        
        !----------------------------------------------------
        ! First loop over each dimension for Lees Edwards b.c.
        ! afterwards search for periodic b.c.
        !----------------------------------------------------
        
        IF ( num_le > 0 ) THEN
           
           DO i = 1, dim
              
              IF ( this%bcdef(2*i-1) == ppm_param_bcdef_LE ) THEN 
                 
                 num = num + 1
                 
                 DO k = 1, dim
                    
                    IF ( k == i ) THEN
                       
                       this%x_image(k,1:num_colloid,num) = &
                            this%x(k,1:num_colloid) - &
                            length(k)
                       
                    ELSE
                       
                       DO j = 1, num_colloid
                          
                          this%x_image(k,j,num) = &
                               MODULO(this%x(k,j) + shear_length(k,2*i-1), &
                               length(k) )
                          
                          this%v_image(k,j,num) = &
                               this%v(k,j,1) + &
                               ( shear_v(k,2*i-1) - shear_v(k,2*i) )
                          
                       END DO ! j = 1, num_colloid
                       
                    END IF ! k == i
                    
                 END DO ! k = 1, dim
                 
              END IF ! bcdef(2*i-1)
              
              IF ( this%bcdef(2*i) == ppm_param_bcdef_LE ) THEN
                 
                 num = num + 1
                 
                 DO k = 1, dim
                    
                    IF ( k == i ) THEN
                       
                       this%x_image(k,1:num_colloid,num) = &
                            this%x(k,1:num_colloid) + &
                            length(k)
                       
                    ELSE
                       
                       DO j = 1, num_colloid
                          
                          this%x_image(k,j,num) = &
                               MODULO(this%x(k,j) + shear_length(k,2*i), &
                               length(k) )
                          
                          this%v_image(k,j,num) = &
                               this%v(k,j,1) + &
                               ( shear_v(k,2*i) - shear_v(k,2*i-1) )
                          
                       END DO ! j = 1, num_colloid
                       
                    END IF ! k == i
                    
                 END DO ! k = 1 , dim
                 
              END IF ! bcdef(2*i)
              
           END DO ! i = 1, dim
           
        END IF ! num_le  > 0
        
        
        !----------------------------------------------------
        ! Loop over each dimension for periodic b.c.
        !----------------------------------------------------
        
        DO i = 1, dim
           
           num_image = num
           
           IF ( this%bcdef(2*i-1) == ppm_param_bcdef_periodic ) THEN
              
              DO j = 1, num_image
                 
                 num = num + 1
                 
                 DO k = 1, dim
                    
                    IF ( k == i ) THEN
                       
                       this%x_image(k,1:num_colloid,num) = &
                            this%x_image(k,1:num_colloid,j) - &
                            length(k)

                    ELSE
                       
                       this%x_image(k,1:num_colloid,num) = &
                            this%x_image(k,1:num_colloid,j)

                    END IF
                    
                 END DO
                 
              END DO ! j = 1, num_image
              
           END IF
           
           IF ( this%bcdef(2*i) == ppm_param_bcdef_periodic ) THEN
              
              DO j = 1, num_image
                 
                 num = num + 1
                 
                 DO k = 1, dim

                    IF ( k == i ) THEN
                       
                       this%x_image(k,1:num_colloid,num) = &
                            this%x_image(k,1:num_colloid,j) + &
                            length(k)
                    
                    ELSE
                       
                       this%x_image(k,1:num_colloid,num) = &
                            this%x_image(k,1:num_colloid,j)
                 
                    END IF
                    
                 END DO
                 
              END DO ! j = 1, num_image
                 
           END IF
              
        END DO ! i = 1, dim
        
        
        num_image = num
        
        IF ( num_image /= this%num_image ) THEN

           PRINT *, "colloid_compute_image : " , &
                "Number of images unexpected : ", &
                num_image
           stat_info = -1
           GOTO 9999
        END IF
        
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Release dynamic memrory
        !----------------------------------------------------
        
        IF (ASSOCIATED(shear_length) ) THEN
           DEALLOCATE(shear_length)
        END IF
        
        IF (ASSOCIATED(shear_v) ) THEN
           DEALLOCATE(shear_v)
        END IF
        
        RETURN
        
      END SUBROUTINE colloid_compute_image
      
