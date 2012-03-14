      SUBROUTINE colloid_initialize_image(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_initialize_image
        !----------------------------------------------------
        !
        ! Purpose     : Initialize the images of all colloids
        !               and save them in the object data
        !               structure for future usage.
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
        INTEGER                         :: num_le,num_peri
        INTEGER                         :: num, num_image
        
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------

        stat_info     = 0
        stat_info_sub = 0
        
        dim         = this%num_dim
        num_colloid = this%num_colloid
        
        !----------------------------------------------------
        ! Check for number of periodic and Lees Edwards b.cs.
        ! The number of image will be related by power of 3.
        !----------------------------------------------------
        
        num_peri = &
             boundary_get_num_peri(this%boundary,stat_info_sub)       
        num_le = &
             boundary_get_num_le(this%boundary,stat_info_sub)
        
        num = num_peri + num_le
        
        num_image = 3**(num/2)
        
        !----------------------------------------------------
        ! Allocate memory for images and initialize.
        !----------------------------------------------------
        
        IF ( ASSOCIATED(this%x_image)) THEN
           DEALLOCATE(this%x_image)
        END IF
        IF ( ASSOCIATED(this%v_image)) THEN
           DEALLOCATE(this%v_image)
        END IF
        
        ALLOCATE(this%x_image(dim,num_colloid,num_image))
        ALLOCATE(this%v_image(dim,num_colloid,num_image))
        
        this%x_image(:,:,:) = 0.0_MK
        this%v_image(:,:,:) = 0.0_MK
        
        this%num_image = num_image
        
        RETURN
        
      END SUBROUTINE colloid_initialize_image
      
