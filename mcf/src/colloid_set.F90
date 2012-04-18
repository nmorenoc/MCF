!--------------------------------------------------
! Subroutine  : colloid_set_*
!--------------------------------------------------
!
! Purpose     : Set routines of Class colloid.
!
! Reference   :
!
! Remark      :
!
! Revisions   : V0.1 01.03.2009, original version.
!
!--------------------------------------------------
! Author      : Xin Bian
! Contact     : xin.bian@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!--------------------------------------------------

      SUBROUTINE colloid_set_tech(this,d_tech,stat_info)
        TYPE(Colloid), INTENT(INOUT)            :: this
        TYPE(Technique),INTENT(IN),TARGET       :: d_tech
        INTEGER, INTENT(OUT)                    :: stat_info

        stat_info = 0
        this%tech => d_tech
        
        RETURN
      END SUBROUTINE colloid_set_tech

      
      SUBROUTINE colloid_set_num_dim(this,d_dim, stat_info)
        !-----------------------------------------
        ! Set the number of dimension of colloids.
        !-----------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_dim
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: num

        num = this%num_colloid
        stat_info = 0
        
        !---------------------------------------
        ! Only 2D, 3D are supported
        !---------------------------------------
        
        IF(d_dim < 2 .OR. d_dim > 3 ) THEN
           PRINT *, "colloid_set_num_dim : ", &
                "Dimension is not supported !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !---------------------------------------
        ! If reset the dimension, 
        ! the memory has to be reallocated
        !---------------------------------------
        IF( d_dim /= this%num_dim ) THEN
           
           this%num_dim = d_dim
           
           IF(ASSOCIATED(this%body_force))THEN
              DEALLOCATE(this%body_force)
           END IF
           ALLOCATE(this%body_force(d_dim))
           
           IF(ASSOCIATED(this%radius))THEN
              DEALLOCATE(this%radius)
           END IF
           ALLOCATE(this%radius(d_dim,num))
           
           IF(ASSOCIATED(this%x))THEN
              DEALLOCATE(this%x)
           END IF
           ALLOCATE(this%x(d_dim,num))
           
           IF(ASSOCIATED(this%v))THEN
              DEALLOCATE(this%v)
           END IF
           ALLOCATE(this%v(d_dim,num,this%integrate_type))
           
#if __DRAG_PART
           IF(ASSOCIATED(this%drag_lub))THEN
              DEALLOCATE(this%drag_lub)
           END IF
           ALLOCATE(this%drag_lub(d_dim,num))

           IF(ASSOCIATED(this%drag_repul))THEN
              DEALLOCATE(this%drag_repul)
           END IF
           ALLOCATE(this%drag_repul(d_dim,num))
#endif
           
           IF(ASSOCIATED(this%drag))THEN
              DEALLOCATE(this%drag)
           END IF
           ALLOCATE(this%drag(d_dim,num))
           
           !-------------------------------------------------
           ! Derived quantities.                               
           !-------------------------------------------------
           
           IF(ASSOCIATED(this%f))THEN
              DEALLOCATE(this%f)
           END IF
           ALLOCATE(this%f(d_dim,num,this%integrate_type))   
           
           IF(ASSOCIATED(this%mom))THEN
              DEALLOCATE(this%mom)
           END IF
           ALLOCATE(this%mom(d_dim,num))
           
           IF(ASSOCIATED(this%mom_tot))THEN
              DEALLOCATE(this%mom_tot)
           END IF
           ALLOCATE(this%mom_tot(d_dim))
           
           !-------------------------------------------------
           ! Physics parameters, boundaries.
           !-------------------------------------------------
           
           IF(ASSOCIATED(this%min_phys)) THEN
              DEALLOCATE(this%min_phys)
           END IF
           ALLOCATE(this%min_phys(d_dim))
           
           IF(ASSOCIATED(this%max_phys)) THEN
              DEALLOCATE(this%max_phys)
           END IF
           ALLOCATE(this%max_phys(d_dim))
           
           IF(ASSOCIATED(this%bcdef)) THEN
              DEALLOCATE(this%bcdef)
           END IF
           ALLOCATE(this%bcdef(2*d_dim))
           
        END IF

9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_num_dim
      

      SUBROUTINE colloid_set_num_colloid(this,d_num,stat_info)
        !-----------------------------------------
        ! Set the number of colloids,
        ! not used
        !-----------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_num
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        
        stat_info = 0
        
        dim = this%num_dim
        
        !---------------------------------------
        ! If reset the number of colloids, 
        ! the memory has to be reallocated
        !---------------------------------------
        
        IF( d_num /= this%num_colloid ) THEN
           
           this%num_colloid = d_num
           
           IF(ASSOCIATED(this%shape))THEN
              DEALLOCATE(this%shape)
           END IF
           ALLOCATE(this%shape(d_num))
           
           IF(ASSOCIATED(this%radius))THEN
              DEALLOCATE(this%radius)
           END IF
           ALLOCATE(this%radius(dim,d_num))
           
           IF(ASSOCIATED(this%freq))THEN
              DEALLOCATE(this%freq)
           END IF
           ALLOCATE(this%freq(d_num))
        
           IF(ASSOCIATED(this%m))THEN
              DEALLOCATE(this%m)
           END IF
           ALLOCATE(this%m(d_num))
           
           IF(ASSOCIATED(this%mmi))THEN
              DEALLOCATE(this%mmi)
           END IF
           ALLOCATE(this%mmi(3,d_num))
           
           IF(ASSOCIATED(this%x))THEN
              DEALLOCATE(this%x)
           END IF
           ALLOCATE(this%x(dim,d_num))
           
           IF(ASSOCIATED(this%v))THEN
              DEALLOCATE(this%v)
           END IF
           ALLOCATE(this%v(dim,d_num,this%integrate_type))
           
#if __DRAG_PART
           IF(ASSOCIATED(this%drag_lub))THEN
              DEALLOCATE(this%drag_lub)
           END IF
           ALLOCATE(this%drag_lub(dim,d_num))
           
           IF(ASSOCIATED(this%drag_repul))THEN
              DEALLOCATE(this%drag_repul)
           END IF
           ALLOCATE(this%drag_repul(dim,d_num))
#endif           
           
           IF(ASSOCIATED(this%drag))THEN
              DEALLOCATE(this%drag)
           END IF
           ALLOCATE(this%drag(dim,d_num))
           
           IF(ASSOCIATED(this%rot_vector))THEN
              DEALLOCATE(this%rot_vector)
           END IF
           ALLOCATE(this%rot_vector(4,d_num))

           IF(ASSOCIATED(this%acc_vector))THEN
              DEALLOCATE(this%acc_vector)
           END IF
           ALLOCATE(this%acc_vector(4,d_num))

           IF(ASSOCIATED(this%rot_matrix))THEN
              DEALLOCATE(this%rot_matrix)
           END IF
           ALLOCATE(this%rot_matrix(3,3,d_num))

           IF(ASSOCIATED(this%acc_matrix))THEN
              DEALLOCATE(this%acc_matrix)
           END IF
           ALLOCATE(this%acc_matrix(3,3,d_num))

           IF(ASSOCIATED(this%theta))THEN
              DEALLOCATE(this%theta)
           END IF           
           ALLOCATE(this%theta(3,d_num))
        
           IF(ASSOCIATED(this%omega))THEN
              DEALLOCATE(this%omega)
           END IF
           
           ALLOCATE(this%omega(3,d_num,this%integrate_type))
           
           IF(ASSOCIATED(this%torque))THEN
              DEALLOCATE(this%torque)
           END IF
           ALLOCATE(this%torque(3,d_num))
           
           IF(ASSOCIATED(this%num_physical_part))THEN
              DEALLOCATE(this%num_physical_part)
           END IF
           ALLOCATE(this%num_physical_part(d_num))
           
           IF(ASSOCIATED(this%num_numerical_part))THEN
              DEALLOCATE(this%num_numerical_part)
           END IF
           ALLOCATE(this%num_numerical_part(d_num))
           
           !-------------------------------------------------
           ! Derived quantities.                               
           !-------------------------------------------------
           
           IF(ASSOCIATED(this%f))THEN
              DEALLOCATE(this%f)
           END IF
           ALLOCATE(this%f(dim,d_num,this%integrate_type))
           
           IF(ASSOCIATED(this%alpha))THEN
              DEALLOCATE(this%alpha)
           END IF
           ALLOCATE(this%alpha(3,d_num,this%integrate_type))
          
           IF(ASSOCIATED(this%k_energy))THEN
              DEALLOCATE(this%k_energy)
           END IF
           ALLOCATE(this%k_energy(d_num))
           
           IF(ASSOCIATED(this%mom))THEN
              DEALLOCATE(this%mom)
           END IF
           ALLOCATE(this%mom(dim,d_num))
           
        END IF
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_num_colloid


      SUBROUTINE colloid_set_adapt_t_coef(this,d_coef,stat_info)
        !-----------------------------------------
        ! Set the coefficient of adaptive time step
        !-----------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_coef
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0

        IF ( d_coef < 0.0_MK  ) THEN
           PRINT *, "colloid_set_adapt_t_coef : ", &
                "Wrong coefficient !"
           stat_info = -1
        END IF
        
        this%adapt_t_coef = d_coef
        
        RETURN
        
      END SUBROUTINE colloid_set_adapt_t_coef


      SUBROUTINE colloid_set_rho(this,d_rho,stat_info)
        !-----------------------------------------
        ! Set the colloid density.
        !-----------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_rho
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0

        IF ( d_rho < 0.0_MK  ) THEN
           PRINT *, "colloid_set_rho : ", &
                "Wrong rho !"
           stat_info = -1
        END IF
        
        this%rho = d_rho
        
        RETURN
        
      END SUBROUTINE colloid_set_rho
      

      SUBROUTINE colloid_set_rho_type(this,d_rho_type, stat_info)
        !-----------------------------------------
        ! Set the colloid density type.
        !-----------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_rho_type
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0

        IF ( d_rho_type < 0 .OR. d_rho_type > 1 ) THEN
           PRINT *, "colloid_set_rho_type : ", &
                "Wrong rho type !"
           stat_info = -1
        END IF
        
        this%rho_type = d_rho_type
        
        RETURN
        
      END SUBROUTINE colloid_set_rho_type

      
      SUBROUTINE colloid_set_translate(this,d_translate,stat_info)
        !---------------------------------------
        ! Set if colloids can translate or not.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        LOGICAL, INTENT(IN)             :: d_translate
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%translate = d_translate
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_translate
      
      
      SUBROUTINE colloid_set_rotate(this,d_rotate,stat_info)
        !---------------------------------------
        ! Set if colloids can rotate or not.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        LOGICAL, INTENT(IN)             :: d_rotate
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%rotate = d_rotate
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_rotate
      

      SUBROUTINE colloid_set_place(this,d_place,stat_info)
        !---------------------------------------
        ! Set the type of noslip condition on
        ! the surface of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_place
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%place = d_place
        
9999    CONTINUE        
        RETURN       
        
      END SUBROUTINE colloid_set_place
      

      SUBROUTINE colloid_set_noslip_type(this,d_noslip_type,stat_info)
        !---------------------------------------
        ! Set the type of noslip condition on
        ! the surface of colloids.
        !---------------------------------------
  
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_noslip_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        this%noslip_type = d_noslip_type
        
9999    CONTINUE        
        RETURN       
        
      END SUBROUTINE colloid_set_noslip_type


      SUBROUTINE colloid_set_body_force_type(this,d_type,stat_info)
        !---------------------------------------
        ! Set the type of body force on
        ! colloids.
        !---------------------------------------
  
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        this%body_force_type = d_type
        
9999    CONTINUE        
        RETURN       
        
      END SUBROUTINE colloid_set_body_force_type


      SUBROUTINE colloid_set_body_force(this,d_body_force,stat_info)
        !---------------------------------------
        ! Set the body force on all colloids
        !---------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: d_body_force
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        
        stat_info = 0
        
        !---------------------------------------
        ! Check if the input dimension matches.
        !---------------------------------------
        
        dim = SIZE(d_body_force,1)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "colloid_set_body_force : ", "Wrong Dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%body_force(1:dim) = d_body_force(1:dim)
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_body_force
      
      
      
      SUBROUTINE colloid_set_cc_lub_type(this,d_type,stat_info)
        !----------------------------------------------------
        ! Set lubrication interaction type between 
        ! colloid and colloid.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF ( d_type < 0 .OR. d_type > 2 ) THEN
           PRINT *, "colloid_set_cc_lub_type : ", &
                "Wrong lub type !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%cc_lub_type = d_type
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cc_lub_type

      
      SUBROUTINE colloid_set_cc_repul_type(this,d_type,stat_info)
        !----------------------------------------------------
        ! Set repulsive force interaction type between 
        ! colloid and colloid.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF ( d_type < 0 .OR. d_type > 2 ) THEN
           PRINT *, "colloid_set_cc_repul_type : ", &
                "Wrong repul type !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%cc_repul_type = d_type
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cc_repul_type


      SUBROUTINE colloid_set_cc_lub_cut_off(this,d_cut_off,stat_info)
        !----------------------------------------------------
        ! Set cut off for lubrication interaction 
        ! between colloid and colloid.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_cut_off
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%cc_lub_cut_off = d_cut_off
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cc_lub_cut_off
      

      SUBROUTINE colloid_set_cc_lub_cut_on(this,d_cut_on,stat_info)
        !----------------------------------------------------
        ! Set cut on for lubrication interaction 
        ! between colloid and colloid.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_cut_on
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%cc_lub_cut_on = d_cut_on
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cc_lub_cut_on

      
       SUBROUTINE colloid_set_cc_repul_cut_off(this,d_cut_off,stat_info)
         !----------------------------------------------------
        ! Set cut off of repulsive interaction 
        ! between colloid and colloid.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_cut_off
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%cc_repul_cut_off = d_cut_off
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cc_repul_cut_off

      
      SUBROUTINE colloid_set_cc_repul_cut_on(this,d_cut_on,stat_info)
        !----------------------------------------------------
        ! Set cut on of repulsive interaction 
        ! between colloid and colloid.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_cut_on
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%cc_repul_cut_on = d_cut_on
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cc_repul_cut_on
      

      SUBROUTINE colloid_set_cc_repul_F0(this,d_F0,stat_info)
        !----------------------------------------------------
        ! Set maximum repulsive force between colloids.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_F0
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%cc_repul_F0 = d_F0
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cc_repul_F0


      SUBROUTINE colloid_set_cw_lub_type(this,d_type, stat_info)
        !----------------------------------------------------
        ! Set lubrication interaction type between 
        ! colloid and wall.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF ( d_type < 0 .OR. d_type > 2 ) THEN
           PRINT *, "colloid_set_cw_lub_type : ", &
                "Wrong lub type !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%cw_lub_type = d_type
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cw_lub_type
      
      
      SUBROUTINE colloid_set_cw_repul_type(this,d_type,stat_info)
        !----------------------------------------------------
        ! Set repulsive force interaction type between 
        ! colloid and wall.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(IN)             :: d_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF ( d_type < 0 .OR. d_type > 2 ) THEN
           PRINT *, "colloid_set_cw_repul_type : ", &
                "Wrong repul type !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%cw_repul_type = d_type
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cw_repul_type
      
      
      SUBROUTINE colloid_set_cw_lub_cut_off(this,d_cut_off, stat_info)
        !----------------------------------------------------
        ! Set cut off for lubrication interaction 
        ! between colloid and wall.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_cut_off
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%cw_lub_cut_off = d_cut_off
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cw_lub_cut_off
      
      
      SUBROUTINE colloid_set_cw_lub_cut_on(this,d_cut_on,stat_info)
        !----------------------------------------------------
        ! Set cut on for lubrication interaction 
        ! between colloid and wall.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_cut_on
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%cw_lub_cut_on = d_cut_on
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cw_lub_cut_on

      
      SUBROUTINE colloid_set_cw_repul_cut_off(this,d_cut_off,stat_info)
        !----------------------------------------------------
        ! Set cut off of repulsive interaction 
        ! between colloid and wall.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_cut_off
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%cw_repul_cut_off = d_cut_off
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cw_repul_cut_off
      
      
      SUBROUTINE colloid_set_cw_repul_cut_on(this,d_cut_on,stat_info)
        !----------------------------------------------------
        ! Set cut on of repulsive interaction 
        ! between colloid and wall.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_cut_on
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%cw_repul_cut_on = d_cut_on
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cw_repul_cut_on


      SUBROUTINE colloid_set_cw_repul_F0(this,d_F0,stat_info)
        !----------------------------------------------------
        ! Set maximum repulsive force between colloid and wall.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_F0
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%cw_repul_F0 = d_F0
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cw_repul_F0

      
      SUBROUTINE colloid_set_shape(this, d_shape,stat_info)
        !---------------------------------------
        ! Set the shape of a colloid object
        !---------------------------------------
  
        TYPE(Colloid), INTENT(INOUT)            :: this
        INTEGER, DIMENSION(:), INTENT(IN)       :: d_shape
        INTEGER, INTENT(OUT)                    :: stat_info
       
        INTEGER                                 :: num

        stat_info = 0
        
        num = SIZE(d_shape,1)
        
        IF (num /= this%num_colloid) THEN
           PRINT *, "colloid_set_shape : ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shape(:) = d_shape(1:num)

9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_shape
      
      
      SUBROUTINE colloid_set_radius(this,d_ra,stat_info)
        !---------------------------------------
        ! Set the radius of a colloid object 
        !---------------------------------------

        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:),INTENT(IN)     :: d_ra
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim,num

        stat_info = 0
        
        dim = SIZE(d_ra,1)
        num = SIZE(d_ra,2)
        
        IF (dim /= this%num_dim) THEN
           PRINT *, "colloid_set_radius : ", &
                "Number of dimension doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF (num /= this%num_colloid) THEN
           PRINT *, "colloid_set_ra : ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%radius(:,:) = d_ra(:,:)

9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_radius
      
      
      SUBROUTINE colloid_set_freq(this,d_freq,stat_info)
        !---------------------------------------
        ! Set frequency of surface roughness.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        INTEGER, DIMENSION(:),INTENT(IN)        :: d_freq
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: num
        
        stat_info = 0
        
        num = SIZE(d_freq,1)
        
        IF ( num /= this%num_colloid ) THEN
           PRINT *, "colloid_set_freq : ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%freq(:) = d_freq(:)
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_freq
      
      
      SUBROUTINE colloid_set_m(this, d_m, stat_info)
        !---------------------------------------
        ! Set the total mass of a colloid object
        !---------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: d_m
        INTEGER, INTENT(OUT)                    :: stat_info

        INTEGER                                 :: num

        
        stat_info = 0
        
        num = SIZE(d_m,1)
  
        IF (num /= this%num_colloid) THEN
           PRINT *, "colloid_set_m : ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%m(:) = d_m(1:num)

9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_m
      
      
      SUBROUTINE colloid_set_mmi(this, d_mmi, stat_info)
        !----------------------------------------------------
        ! Set the mass momentum inertia of colloids
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_mmi
        INTEGER, INTENT(OUT)                    :: stat_info

        INTEGER                                 :: dim,num

        
        stat_info = 0
        
        dim = SIZE(d_mmi,1)
        num = SIZE(d_mmi,2)
        
        IF (dim /= 3) THEN
           PRINT *, "colloid_set_mmi : ", &
                "Number of dimensions doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF (num /= this%num_colloid) THEN
           PRINT *, "colloid_set_mmi: ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%mmi(:,:) = d_mmi(:,:)

9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_mmi


      SUBROUTINE colloid_set_x(this,d_x,stat_info)
        !----------------------------------------------------
        ! Set the positions of colloids.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_x
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        INTEGER                                 :: num

        stat_info = 0
                
        !----------------------------------------------------
        ! Check if the input position's dimension
        ! and num of colloids match.
        !---------------------------------------------------
        dim = SIZE(d_x,1)
        num = SIZE(d_x,2)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "colloid_set_x: ", "Wrong Dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( num /= this%num_colloid) THEN
           PRINT *, "colloid_set_x: ", "Wrong number of colloids !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%x(:,:) = d_x(1:dim,1:num)

9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_x
      

      SUBROUTINE colloid_set_v(this,d_v,stat_info)
        !----------------------------------------------------
        ! Set the velocity of colloids.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK),DIMENSION(:,:,:), INTENT(IN)   :: d_v
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        INTEGER                                 :: num
        INTEGER                                 :: itype

        stat_info = 0
        
        !----------------------------------------------------
        ! Check if the input velocity's dimension
        ! and num of colloids match.
        !----------------------------------------------------
        
        dim   = SIZE(d_v,1)
        num   = SIZE(d_v,2)
        itype = SIZE(d_v,3)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "colloid_set_v: ", "Wrong Dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( num /= this%num_colloid) THEN
           PRINT *, "colloid_set_v: ", "Wrong number of colloids !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( itype /= this%integrate_type) THEN
           PRINT *, "colloid_set_v: ", "Wrong integration type!"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%v(1:dim,1:num,1:itype) = d_v(1:dim,1:num,1:itype)
        
9999    CONTINUE
        RETURN       
        
      END SUBROUTINE colloid_set_v
      
      
#if __DRAG_PART
      SUBROUTINE colloid_set_drag_lub(this,d_drag,stat_info)
        !----------------------------------------------------
        ! Set the drag/force lubrication part
        ! exerted on a colloid object
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        INTEGER                                 :: num
        
        stat_info = 0
        
        !----------------------------------------------------
        ! Check if the input force's dimension and
        ! num of colloids match.
        !----------------------------------------------------
        
        dim = SIZE(d_drag,1)
        num = SIZE(d_drag,2)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "colloid_set_drag_lub: ", "Wrong Dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( num /= this%num_colloid) THEN
           PRINT *, "colloid_set_drag_lub: ", "Wrong number of colloids !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%drag_lub(:,:) = d_drag(1:dim,1:num)

9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_drag_lub

      
      SUBROUTINE colloid_set_drag_repul(this,d_drag,stat_info)
        !----------------------------------------------------
        ! Set the drag/force repulsive part
        ! exerted on a colloid object
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        INTEGER                                 :: num
        
        
        stat_info = 0
        
        !---------------------------------------
        ! Check if the input force's dimension and
        ! num of colloids match.
        !---------------------------------------
        
        dim = SIZE(d_drag,1)
        num = SIZE(d_drag,2)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "colloid_set_drag: ", "Wrong Dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( num /= this%num_colloid) THEN
           PRINT *, "colloid_set_drag: ", "Wrong number of colloids !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%drag_repul(:,:) = d_drag(1:dim,1:num)

9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_drag_repul
#endif

      
      SUBROUTINE colloid_set_drag(this,d_drag,stat_info)
        !----------------------------------------------------
        ! Set the drag/force exerted on a colloid object
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        INTEGER                                 :: num
        
        
        stat_info = 0
        
        !----------------------------------------------------
        ! Check if the input force's dimension and
        ! num of colloids match.
        !----------------------------------------------------
        
        dim = SIZE(d_drag,1)
        num = SIZE(d_drag,2)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "colloid_set_drag: ", "Wrong Dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( num /= this%num_colloid) THEN
           PRINT *, "colloid_set_drag: ", "Wrong number of colloids !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%drag(:,:) = d_drag(1:dim,1:num)

9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_drag
      
      
      SUBROUTINE colloid_add_drag(this, d_drag,stat_info)
        !----------------------------------------------------
        ! Accumulate the drag/force on colloids.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:),INTENT(IN)     :: d_drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        INTEGER                                 :: num

        stat_info = 0
        
        dim = SIZE(d_drag,1)
        num = SIZE(d_drag,2)
        
        !----------------------------------------------------
        ! check if the input force's dimension matches
        !----------------------------------------------------
        
        IF (dim /= this%num_dim) THEN
           PRINT *, "colloid_add_drag: ", &
                "Dimension doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        IF (num /= this%num_colloid) THEN
           PRINT *, "colloid_add_drag: ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%drag(1:dim,1:num) =  &
             this%drag(1:dim,1:num) + d_drag(1:dim,1:num)
        
9999    CONTINUE
        RETURN       
        
      END SUBROUTINE colloid_add_drag


      SUBROUTINE colloid_set_rotation_vector(this,d_vector,stat_info)
        !----------------------------------------------------
        ! Set rotation vector
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:),INTENT(IN)     :: d_vector
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: num
        
        stat_info = 0
        
        num = SIZE(d_vector,1)
        
        IF (num /= 4) THEN
           PRINT *, "colloid_set_rotation_vector: ", &
                "Number of dimension doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF

        num = SIZE(d_vector,2)
        
        IF (num /= this%num_colloid) THEN
           PRINT *, "colloid_set_rotation_vector: ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%rot_vector(:,:) = d_vector(:,:)
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_rotation_vector
      
            
      SUBROUTINE colloid_set_accumulation_vector(this,d_vector,stat_info)
        !----------------------------------------------------
        ! Set accumulative rotation vector.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:),INTENT(IN)     :: d_vector
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num,i
        REAL(MK)                                :: len
        
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------
        
        stat_info = 0
        stat_info_sub = 0
        
        !----------------------------------------------------
        ! check dimension
        !----------------------------------------------------
        
        num = SIZE(d_vector,1)
        
        IF (num /= 4) THEN
           PRINT *, "colloid_set_accumulation_vector : ", &
                "Number of dimension doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        num = SIZE(d_vector,2)
        
        IF (num /= this%num_colloid) THEN
           PRINT *, "colloid_set_acccumulation_vector : ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Normalize rotation axis.
        ! If axis is a zero vector, set rotation angle=0.
        !----------------------------------------------------
        
        DO i = 1, this%num_colloid
           
           len = SQRT(DOT_PRODUCT(d_vector(1:3,i), d_vector(1:3,i)))
           
           IF ( len < mcf_machine_zero ) THEN
              
              this%acc_vector(1,i) = 1.0_MK
              this%acc_vector(2,i) = 0.0_MK
              this%acc_vector(3,i) = 0.0_MK
              this%acc_vector(4,i) = 0.0_MK
              
           ELSE

              this%acc_vector(1:3,i) = d_vector(1:3,i) / len
              this%acc_vector(4,i)   = d_vector(4,i)
              
           END IF
           
        END DO        
        
        
        !----------------------------------------------------
        ! Use initial accumulation vector to initialize
        ! accumulation matrix.
        !----------------------------------------------------
        
        CALL colloid_init_accumulation_matrix(this,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "colloid_set_accumulation_vector: ", &
                "compute accumulation matrix failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_accumulation_vector
      

      SUBROUTINE colloid_set_rotation_matrix(this,d_matrix,stat_info)
        !----------------------------------------------------
        ! Set current rotation matrix.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:,:),INTENT(IN)   :: d_matrix
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim1,dim2, num
        
        stat_info = 0
        
        dim1 = SIZE(d_matrix,1)
        dim2 = SIZE(d_matrix,2)
        num  = SIZE(d_matrix,3)
        
        IF ( dim1 /= 3 .OR. &
             dim2 /= 3 ) THEN
           PRINT *, "colloid_set_rotation_matrix: ", &
                "Number of dimension doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF

        IF ( num /= this%num_colloid ) THEN
           PRINT *, "colloid_set_rotation_matrix: ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%rot_matrix(:,:,:) = d_matrix(:,:,:)
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_rotation_matrix
      

      SUBROUTINE colloid_set_accumulation_matrix(this,d_matrix,stat_info)
        !----------------------------------------------------
        ! Set accumulative rotation matrix.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:,:),INTENT(IN)   :: d_matrix
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim1,dim2, num
        
        stat_info = 0
        
        dim1 = SIZE(d_matrix,1)
        dim2 = SIZE(d_matrix,2)
        num  = SIZE(d_matrix,3)
        
        IF ( dim1 /= 3 .OR. &
             dim2 /= 3 ) THEN
           PRINT *, "colloid_set_accumulation_matrix : ", &
                "Number of dimension doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF

        IF ( num /= this%num_colloid ) THEN
           PRINT *, "colloid_set_accumulation_matrix : ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%acc_matrix(:,:,:) = d_matrix(:,:,:)
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_accumulation_matrix
      
      
      SUBROUTINE colloid_set_theta(this,d_theta,stat_info)
        !----------------------------------------------------
        ! Set rotated angle.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:),INTENT(IN)       :: d_theta
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: num
        
        stat_info = 0
        
        num = SIZE(d_theta,2)
        
        IF (num /= this%num_colloid) THEN
           PRINT *, "colloid_set_theta : ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%theta(:,:) = d_theta(:,:)
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE colloid_set_theta


      SUBROUTINE colloid_set_omega(this,d_omega,stat_info)
        !----------------------------------------------------
        ! Set the rotating velocity of colloids.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:,:), INTENT(IN)  :: d_omega
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        INTEGER                                 :: num
        INTEGER                                 :: itype
        
        stat_info = 0
        
        !----------------------------------------------------
        ! Check if the input rotating velocity's
        ! dimension is 3.
        !----------------------------------------------------
        dim = SIZE(d_omega,1)
        num = SIZE(d_omega,2)
        itype = SIZE(d_omega,3)
        
        IF( dim /= 3) THEN
           PRINT *, "colloid_set_oemga : ", "Wrong Dimension!"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( num /= this%num_colloid) THEN
           PRINT *, "colloid_set_omega : ", "Wrong number of colloids!"
           stat_info = -1
           GOTO 9999
        END IF

        IF( itype /= this%integrate_type) THEN
           PRINT *, "colloid_set_omega : ", "Wrong integration type!"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%omega(1:3,1:num,1:itype) = d_omega(1:3,1:num,1:itype)
        
9999    CONTINUE
        RETURN       
        
      END SUBROUTINE colloid_set_omega
      
      
      SUBROUTINE colloid_set_torque(this,d_torque,stat_info)
        !----------------------------------------------------
        ! Set the torque of colloids.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: d_torque
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        INTEGER                                 :: num

        stat_info = 0
        
        !----------------------------------------------------
        ! Check if the input torque's
        ! dimension is 3.
        !----------------------------------------------------
        dim = SIZE(d_torque,1)
        num = SIZE(d_torque,2)
        
        IF( dim /= 3) THEN
           PRINT *, "colloid_set_torque : ", "Wrong Dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( num /= this%num_colloid) THEN
           PRINT *, "colloid_set_torque : ", "Wrong number of colloids !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%torque(1:3,1:num) = d_torque(1:3,1:num)
        
9999    CONTINUE
        RETURN       
        
      END SUBROUTINE colloid_set_torque
      

      SUBROUTINE colloid_set_num_physical_part(this,&
           d_num_physical_part,stat_info)
        !----------------------------------------------------
        ! Set the number of physical particles,
        ! which consitute colloids.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        INTEGER, DIMENSION(:), INTENT(IN)       :: d_num_physical_part
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: num, i
        

        stat_info = 0
        
        num = SIZE(d_num_physical_part,1)
  
        IF (num /= this%num_colloid) THEN
           PRINT *, "colloid_set_num_physical_part : ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%num_physical_part(:) = &
             d_num_physical_part(1:num)
        
        this%num_physical_part_tot = 0
        
        DO i =1, this%num_colloid
           
           this%num_physical_part_tot =  &
                this%num_physical_part_tot + &
                this%num_physical_part(i)                   
        END DO
        
9999    CONTINUE  
        
        RETURN       
        
      END SUBROUTINE colloid_set_num_physical_part
      

      SUBROUTINE colloid_set_num_numerical_part(this,&
           d_num_numerical_part,stat_info)
        !----------------------------------------------------
        ! Set the number of numerical particles,
        ! which consitute colloids.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        INTEGER, DIMENSION(:), INTENT(IN)       :: d_num_numerical_part
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: num,i

        stat_info = 0
        
        num = SIZE(d_num_numerical_part,1)
  
        IF (num /= this%num_colloid) THEN
           PRINT *, "colloid_set_num_numerical_part : ", &
                "Number of colloids doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        this%num_numerical_part(:) = &
             d_num_numerical_part(:)
        
        this%num_numerical_part_tot = 0
        
        DO i =1, this%num_colloid
           
           this%num_numerical_part_tot =  &
                this%num_numerical_part_tot + &
                this%num_numerical_part(i)                   
        END DO
        
9999    CONTINUE        
        RETURN       
        
      END SUBROUTINE colloid_set_num_numerical_part
      
      
      SUBROUTINE colloid_set_min_phys(this,d_min_phys,stat_info)
        
        TYPE(Colloid), INTENT(INOUT)             :: this
        REAL(MK), DIMENSION(:)                   :: d_min_phys
        INTEGER, INTENT(OUT)                     :: stat_info
        
        INTEGER                                  :: dim
        
        stat_info = 0
        dim = SIZE(d_min_phys)
        IF (dim /= this%num_dim ) THEN
           PRINT *, "colloid_set_min_phys : ", "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%min_phys(1:dim) = d_min_phys(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_min_phys

      
      SUBROUTINE colloid_set_max_phys(this,d_max_phys,stat_info)
        
        TYPE(Colloid), INTENT(INOUT)             :: this
        REAL(MK), DIMENSION(:)                   :: d_max_phys
        INTEGER, INTENT(OUT)                     :: stat_info
        
        INTEGER                                  :: dim
        
        stat_info = 0
        dim = SIZE(d_max_phys)
        IF (dim /= this%num_dim ) THEN
           PRINT *, "colloid_set_max_phys : ", "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%max_phys(1:dim) = d_max_phys(1:dim)

9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_max_phys
      

      SUBROUTINE colloid_set_min_phys_t(this,d_min_phys_t,stat_info)
        
        TYPE(Colloid), INTENT(INOUT)             :: this
        REAL(MK), DIMENSION(:)                   :: d_min_phys_t
        INTEGER, INTENT(OUT)                     :: stat_info
        
        INTEGER                                  :: dim
        
        stat_info = 0
        
        dim = SIZE(d_min_phys_t)
        
        IF (dim /= this%num_dim ) THEN
           PRINT *, "colloid_set_min_phys_t : ", "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%min_phys_t(1:dim) = d_min_phys_t(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_min_phys_t

      
      SUBROUTINE colloid_set_max_phys_t(this,d_max_phys_t,stat_info)
        
        TYPE(Colloid), INTENT(INOUT)             :: this
        REAL(MK), DIMENSION(:)                   :: d_max_phys_t
        INTEGER, INTENT(OUT)                     :: stat_info
        
        INTEGER                                  :: dim
        
        stat_info = 0

        dim = SIZE(d_max_phys_t)
        
        IF (dim /= this%num_dim ) THEN
           PRINT *, "colloid_set_max_phys_t : ", "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%max_phys_t(1:dim) = d_max_phys_t(1:dim)

9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_max_phys_t
      
       
      SUBROUTINE colloid_set_cut_off(this,d_cut_off,stat_info)
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_cut_off
        INTEGER, INTENT(OUT)            :: stat_info
        
        REAL(MK)                        :: din
        INTEGER                         :: stat_info_sub
        
        stat_info     = 0
        stat_info_sub = 0
        
        IF (d_cut_off <= 0.0_MK) THEN
           PRINT *, "colloid_set_cut_off : ", &
                "cut off should be non-negative !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%cut_off = d_cut_off
        
        din = (1.0_MK+mcf_colloid_in_layer_coeff)*d_cut_off

        IF ( this%rho_type == mcf_colloid_rho_type_dynamic ) THEN
           
           din = din * 2.0_MK
           
        END IF

        CALL colloid_set_din(this,din, stat_info_sub)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_cut_off
      
      
      SUBROUTINE colloid_set_dout(this,d_dout,stat_info)
        !---------------------------------------
        ! Set the minimal distance of a fluid
        ! particle from the surface of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_dout
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        IF (d_dout < 0.0_MK) THEN
           PRINT *, "colloid_set_dout : ", &
                "dout should be positive !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%dout = d_dout
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_dout
      
      
      SUBROUTINE colloid_set_din(this,d_din, stat_info)
        !---------------------------------------
        ! Set the maximum distance of a boundary
        ! particle from the surface of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_din
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        IF ( d_din < 0.0_MK ) THEN
           PRINT *, "colloid_set_din : ", &
                "din should be positive !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%din = d_din
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_din
      

      SUBROUTINE colloid_set_eta(this,d_eta,stat_info)
        !-----------------------------------------
        ! Set the fluid dynamics viscosity.
        !-----------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        REAL(MK), INTENT(IN)            :: d_eta
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0

        IF ( d_eta < 0.0_MK ) THEN
           PRINT *, "colloid_set_eta : ", &
                "Wrong eta !"
           stat_info = -1
        END IF
        
        this%eta = d_eta
        
        RETURN
        
      END SUBROUTINE colloid_set_eta

      
      SUBROUTINE colloid_set_bcdef(this,d_bcdef,stat_info)

        TYPE(Colloid), INTENT(INOUT)          :: this
        INTEGER, DIMENSION(:)                 :: d_bcdef
        INTEGER, INTENT(OUT)                  :: stat_info
        
        
        INTEGER                               :: bcdef_dim
        
        stat_info = 0
        
        bcdef_dim = SIZE(d_bcdef)
        
        IF( bcdef_dim /= 2*this%num_dim ) THEN
           PRINT *, "colloid_set_bcdef : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
           
        END IF
        
        this%bcdef(1:bcdef_dim) = d_bcdef(1:bcdef_dim)
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  colloid_set_bcdef

      
      SUBROUTINE colloid_set_boundary(this,d_boundary,stat_info)
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        TYPE(Boundary), TARGET          :: d_boundary
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_dim
        
        stat_info     = 0
        stat_info_sub = 0
        
        num_dim = &
             boundary_get_num_dim(d_boundary,stat_info_sub)
        
        IF(num_dim /= this%num_dim) THEN
           PRINT *, "colloid_set_boundary : ",&
                "Boundarys dimension doesn't match with colloid !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%boundary => d_boundary
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_set_boundary
      

#if 0
      SUBROUTINE colloid_set_image(this,stat_info)
        !----------------------------------------------------
        ! Remark      :
        !               In case of periodic or Lees-Edwards
        !               boundaries, the images of colloid's
        !               center has to be taken into account
        !               to decide if a potential particle 
        !               is inside the geometry of a colloid.
        !
        !               For 2D, maximum 3**2=9=8 images;
        !               For 3D, maximum 3**3=27=26 images
        !               (including itself).
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)          :: this
        INTEGER, INTENT(OUT)                  :: stat_info
        
        
        INTEGER                               :: num_dim
        INTEGER                               :: num
        INTEGER, DIMENSION(3)                 :: istart
        INTEGER, DIMENSION(3)                 :: iend
        INTEGER                               :: i,j,k       
        REAL(MK), DIMENSION(3,27)             :: image
        
        
        
        stat_info = 0

        
        num_dim = this%num_dim
        
        istart(:) = 0
        iend(:)   = 0
        
        DO i = 1, num_dim
           
           !-------------------------------------------------
           ! If one side is periodic, 
           ! the opposite must be periodic also.
           ! If one side is Lees-Edwards, 
           ! the opposite must be Lees-Edwards also.
           !-------------------------------------------------
           
           IF ( this%bcdef(2*i-1) == ppm_param_bcdef_periodic .OR. &
                this%bcdef(2*i-1) == ppm_param_bcdef_LE ) THEN
              
              istart(i) = -1
              iend(i)   = 1
              
           END IF
           
        END DO
        
        num = 1
        image(1:num_dim,1) = 0.0_MK
        
        DO k = istart(3), iend(3)
           
           DO j = istart(2), iend(2)
              
              DO i = istart(1), iend(1)
                 
                 !-------------------------------------------
                 ! Box(cell) itself doesn't count,
                 ! we exclude.
                 !-------------------------------------------
                 
                 IF ( i == 0 .AND. &
                      j == 0 .AND. &
                      k == 0 ) THEN
                    
                    CYCLE
                    
                 END IF
                 
                 num = num + 1
                 image(1, num) = i * &
                      ( this%max_phys(1) -this%min_phys(1) )
                 image(2, num) = j * &
                      ( this%max_phys(2) -this%min_phys(2) )
                 
                 IF ( num_dim == 3 ) THEN
                    image(3, num) = k * &
                         ( this%max_phys(3) -this%min_phys(3) )
                 END IF
                 
              END DO ! i
              
           END DO ! j
           
        END DO !  k
        
        this%num_image = num
        
        IF (ASSOCIATED(this%image)) THEN
           DEALLOCATE(this%image)
        END IF
        
        ALLOCATE(this%image(num_dim,num))
        
        this%image(1:num_dim,1:num) = &
             image(1:num_dim,1:num)
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  colloid_set_image
#endif
