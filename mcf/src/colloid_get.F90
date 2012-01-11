!--------------------------------------------------
! Subroutine  : colloid_get_*
!--------------------------------------------------
!
! Purpose     : Get routines of Class Colloid.
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
!-------------------------------------------------

      INTEGER FUNCTION colloid_get_num_dim(this, stat_info)
        !---------------------------------------
        !   Return the num of dimension.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_num_dim = this%num_dim
        
        RETURN

      END FUNCTION colloid_get_num_dim


      INTEGER FUNCTION colloid_get_num_colloid(this, stat_info)
        !---------------------------------------
        !   Return the num of colloids.
        !---------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_num_colloid = this%num_colloid
        
        RETURN
        
      END FUNCTION colloid_get_num_colloid
      
      
      REAL(MK) FUNCTION colloid_get_rho(this, stat_info)
        !---------------------------------------
        ! Return the density of colloids.
        !---------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_rho = this%rho
        
        RETURN
        
      END FUNCTION colloid_get_rho


      INTEGER FUNCTION colloid_get_rho_type(this, stat_info)
        !---------------------------------------
        ! Return the density type of colloids.
        !---------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_rho_type = this%rho_type
        
        RETURN
        
      END FUNCTION colloid_get_rho_type

      
      LOGICAL FUNCTION colloid_get_translate(this,stat_info)
        !--------------------------------------------
        ! Return if the colloids can translate or not.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        colloid_get_translate = this%translate
        
        RETURN
        
      END FUNCTION  colloid_get_translate

      
      LOGICAL FUNCTION colloid_get_rotate(this,stat_info)
        !--------------------------------------------
        ! Return if the colloids can rotate or not.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        colloid_get_rotate = this%rotate
        
        RETURN
        
      END FUNCTION colloid_get_rotate

      
      INTEGER FUNCTION colloid_get_place(this,stat_info)
        !--------------------------------------------
        ! Return how to place boundary particles.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        colloid_get_place = this%place
        
        RETURN
        
      END FUNCTION colloid_get_place
      
      
      INTEGER FUNCTION colloid_get_noslip_type(this,stat_info)
        !--------------------------------------------
        ! Return how to model noslip boundary.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        colloid_get_noslip_type = this%noslip_type
        
        RETURN
        
      END FUNCTION colloid_get_noslip_type


      INTEGER FUNCTION colloid_get_body_force_type(this,stat_info)
        !--------------------------------------------
        ! Return body force type of colloid
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        colloid_get_body_force_type = this%body_force_type
        
        RETURN
        
      END FUNCTION colloid_get_body_force_type

      
      SUBROUTINE colloid_get_body_force(this,d_body_force,stat_info)
        !---------------------------------------
        ! Return the body force of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), POINTER         :: d_body_force
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_body_force)) THEN
           DEALLOCATE(d_body_force)
        END IF
        
        ALLOCATE(d_body_force(this%num_dim))
        
        d_body_force(:) = &
             this%body_force(1:this%num_dim)
        
        RETURN       

      END SUBROUTINE colloid_get_body_force
      
      
      INTEGER FUNCTION colloid_get_cc_lub_type(this, stat_info)
        !----------------------------------------------------
        ! Return lubrication interaction type of colloid-colloid.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cc_lub_type = this%cc_lub_type
        
        RETURN
        
      END FUNCTION colloid_get_cc_lub_type


      INTEGER FUNCTION colloid_get_cc_repul_type(this, stat_info)
        !----------------------------------------------------
        ! Return repulsion interaction type of colloid-colloid.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cc_repul_type = this%cc_repul_type
        
        RETURN
        
      END FUNCTION colloid_get_cc_repul_type
      
      
      REAL(MK) FUNCTION colloid_get_cc_lub_cut_off(this, stat_info)
        !----------------------------------------------------
        ! Return cut off for lubrication interaction
        ! between colloid-colloid.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cc_lub_cut_off = &
             this%cc_lub_cut_off
        
        RETURN
        
      END FUNCTION colloid_get_cc_lub_cut_off


      REAL(MK) FUNCTION colloid_get_cc_lub_cut_on(this, stat_info)
        !----------------------------------------------------
        ! Return cut on for lubrication interaction
        ! between colloid-colloid.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cc_lub_cut_on = &
             this%cc_lub_cut_on
        
        RETURN
        
      END FUNCTION colloid_get_cc_lub_cut_on
      
      
      REAL(MK) FUNCTION colloid_get_cc_repul_cut_off(this, stat_info)
        !----------------------------------------------------
        ! Return cut off for repulsive interaction
        ! between colloid-colloid.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cc_repul_cut_off = &
             this%cc_repul_cut_off
        
        RETURN
        
      END FUNCTION colloid_get_cc_repul_cut_off

      
      REAL(MK) FUNCTION colloid_get_cc_repul_cut_on(this, stat_info)
        !----------------------------------------------------
        ! Return cut on for repulsive interaction
        ! between colloid-colloid.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cc_repul_cut_on = &
             this%cc_repul_cut_on
        
        RETURN
        
      END FUNCTION colloid_get_cc_repul_cut_on


      REAL(MK) FUNCTION colloid_get_cc_repul_F0(this, stat_info)
        !----------------------------------------------------
        ! Return maximum repulsive force
        ! between colloid-colloid.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cc_repul_F0 = &
             this%cc_repul_F0
        
        RETURN
        
      END FUNCTION colloid_get_cc_repul_F0

      
      INTEGER FUNCTION colloid_get_cw_lub_type(this, stat_info)
        !----------------------------------------------------
        ! Return lubrication interaction type of colloid-wall.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cw_lub_type = this%cw_lub_type
        
        RETURN
        
      END FUNCTION colloid_get_cw_lub_type

      
      INTEGER FUNCTION colloid_get_cw_repul_type(this, stat_info)
        !----------------------------------------------------
        ! Return repulsion interaction type of colloid-wall.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cw_repul_type = this%cw_repul_type
        
        RETURN
        
      END FUNCTION colloid_get_cw_repul_type


      REAL(MK) FUNCTION colloid_get_cw_lub_cut_off(this, stat_info)
        !----------------------------------------------------
        ! Return cut off for lubrication interaction
        ! between colloid-wall.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cw_lub_cut_off = &
             this%cw_lub_cut_off
        
        RETURN
        
      END FUNCTION colloid_get_cw_lub_cut_off

      
      REAL(MK) FUNCTION colloid_get_cw_lub_cut_on(this, stat_info)
        !----------------------------------------------------
        ! Return cut on for lubrication interaction
        ! between colloid-wall.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cw_lub_cut_on = &
             this%cw_lub_cut_on
        
        RETURN
        
      END FUNCTION colloid_get_cw_lub_cut_on
      
      
      REAL(MK) FUNCTION colloid_get_cw_repul_cut_off(this, stat_info)
        !----------------------------------------------------
        ! Return cut off for repulsive interaction
        ! between colloid-wall.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cw_repul_cut_off = &
             this%cw_repul_cut_off
        
        RETURN
        
      END FUNCTION colloid_get_cw_repul_cut_off

      
      REAL(MK) FUNCTION colloid_get_cw_repul_cut_on(this, stat_info)
        !----------------------------------------------------
        ! Return cut on for repulsive interaction
        ! between colloid-wall.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cw_repul_cut_on = &
             this%cw_repul_cut_on
        
        RETURN
        
      END FUNCTION colloid_get_cw_repul_cut_on
      

      REAL(MK) FUNCTION colloid_get_cw_repul_F0(this, stat_info)
        !----------------------------------------------------
        ! Return maximum repulsive force
        ! between colloid-wall.
        !----------------------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_cw_repul_F0 = &
             this%cw_repul_F0
        
        RETURN
        
      END FUNCTION colloid_get_cw_repul_F0

      
      SUBROUTINE colloid_get_shape(this, d_shape,stat_info)
        !---------------------------------------
        ! Return the shape of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, DIMENSION(:), POINTER  :: d_shape        
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_shape))THEN
           DEALLOCATE(d_shape)
        END IF

        ALLOCATE(d_shape(this%num_colloid))
        
        d_shape(:) = this%shape(1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_shape
      
      
      SUBROUTINE colloid_get_radius(this,d_ra,stat_info)
        !---------------------------------------
        !  Return the radius of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_ra
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_ra))THEN
           DEALLOCATE(d_ra)
        END IF
        
        ALLOCATE(d_ra(this%num_dim,this%num_colloid))
        
        d_ra(:,:) = this%radius(:,:)
        
        RETURN       
        
      END SUBROUTINE colloid_get_radius

      
      SUBROUTINE colloid_get_freq(this,d_freq,stat_info)
        !---------------------------------------
        ! Return the frequency of surface roughness.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, DIMENSION(:), POINTER  :: d_freq
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_freq))THEN
           DEALLOCATE(d_freq)
        END IF

        ALLOCATE(d_freq(this%num_colloid))
        
        d_freq(:) = this%freq(1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_freq
      
     
      SUBROUTINE colloid_get_m(this,d_m, stat_info)
        !---------------------------------------
        ! Return the total mass of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        REAL(MK), DIMENSION(:), POINTER :: d_m
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_m))THEN
           DEALLOCATE(d_m)
        END IF
        
        ALLOCATE(d_m(this%num_colloid))
        
        d_m(:) = this%m(1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_m


      SUBROUTINE colloid_get_mmi(this,d_mmi, stat_info)
        !---------------------------------------
        ! Return mass momentum inertia of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        REAL(MK),DIMENSION(:,:),POINTER :: d_mmi
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_mmi))THEN
           DEALLOCATE(d_mmi)
        END IF
        
        ALLOCATE(d_mmi(3,this%num_colloid))
        
        d_mmi(:,:) = this%mmi(1:3,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_mmi


      SUBROUTINE colloid_get_x(this,d_x,stat_info)
        !---------------------------------------
        ! Return the positions of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_x
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0

        IF(ASSOCIATED(d_x)) THEN
           DEALLOCATE(d_x)
        END IF
        
        ALLOCATE(d_x(this%num_dim,this%num_colloid))
        
        d_x(:,:) = &
             this%x(1:this%num_dim,1:this%num_colloid)
        
        RETURN
        
      END SUBROUTINE colloid_get_x
      
      
      SUBROUTINE colloid_get_v(this,d_v,stat_info)
        !---------------------------------------
        ! Return the velocity of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_v
        INTEGER, INTENT(OUT)                    :: stat_info
        
        
        stat_info = 0
        
        IF(ASSOCIATED(d_v)) THEN           
           DEALLOCATE(d_v)
        END IF
        
        ALLOCATE(d_v(this%num_dim,this%num_colloid))
        
        d_v(:,:) = this%v(1:this%num_dim,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_v
      
      
#if __DRAG_PART
      SUBROUTINE colloid_get_drag_lub(this,d_drag,stat_info)
        !--------------------------------------------
        ! Return the drag/force lubrication part
        ! exerted on colloids.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_drag)) THEN           
           DEALLOCATE(d_drag)
        END IF
        
        ALLOCATE(d_drag(this%num_dim,this%num_colloid))
        
        d_drag(:,:) = &
             this%drag_lub(1:this%num_dim,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_drag_lub


      SUBROUTINE colloid_get_drag_repul(this,d_drag,stat_info)
        !--------------------------------------------
        ! Return the drag/force repulsive/disperse
        ! part exerted on colloids.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_drag)) THEN           
           DEALLOCATE(d_drag)
        END IF

        ALLOCATE(d_drag(this%num_dim,this%num_colloid))
        
        d_drag(:,:) = &
             this%drag_repul(1:this%num_dim,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_drag_repul
#endif
      
      
      SUBROUTINE colloid_get_drag(this,d_drag,stat_info)
        !--------------------------------------------
        ! Return the drag/force exerted on colloids.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_drag)) THEN           
           DEALLOCATE(d_drag)
        END IF

        ALLOCATE(d_drag(this%num_dim,this%num_colloid))
        
        d_drag(:,:) = &
             this%drag(1:this%num_dim,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_drag


      SUBROUTINE colloid_get_rotation_vector(this,d_vector,stat_info)
        !---------------------------------------
        ! Return rotation vector.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)         :: this
        REAL(MK), DIMENSION(:,:), POINTER :: d_vector
        INTEGER, INTENT(OUT)              :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_vector))THEN
           DEALLOCATE(d_vector)
        END IF
        
        ALLOCATE(d_vector(4,this%num_colloid))
        
        d_vector(:,:) = this%rot_vector(:,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_rotation_vector
      

      SUBROUTINE colloid_get_accumulation_vector(this,&
           d_vector,stat_info)
        !---------------------------------------
        ! Return accumulative rotation vector.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)         :: this
        REAL(MK), DIMENSION(:,:), POINTER :: d_vector
        INTEGER, INTENT(OUT)              :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_vector))THEN
           DEALLOCATE(d_vector)
        END IF
        
        ALLOCATE(d_vector(4,this%num_colloid))
        
        d_vector(:,:) = this%acc_vector(:,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_accumulation_vector

      
      SUBROUTINE colloid_get_rotation_matrix(this,d_matrix,stat_info)
        !---------------------------------------
        ! Return current roation matrix
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK),DIMENSION(:,:,:),POINTER       :: d_matrix
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_matrix))THEN
           DEALLOCATE(d_matrix)
        END IF
        
        ALLOCATE(d_matrix(3,3,this%num_colloid))
        
        d_matrix(:,:,:) = this%rot_matrix(:,:,:)
        
        RETURN       
        
      END SUBROUTINE colloid_get_rotation_matrix


      SUBROUTINE colloid_get_accumulation_matrix(this,d_matrix,stat_info)
        !---------------------------------------
        ! Return accumulative rotation matrix.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK),DIMENSION(:,:,:),POINTER       :: d_matrix
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_matrix))THEN
           DEALLOCATE(d_matrix)
        END IF
        
        ALLOCATE(d_matrix(3,3,this%num_colloid))
        
        d_matrix(:,:,:) = this%acc_matrix(:,:,:)
        
        RETURN       
        
      END SUBROUTINE colloid_get_accumulation_matrix


      SUBROUTINE colloid_get_theta(this,d_theta,stat_info)
        !---------------------------------------
        ! Return roated angle.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)         :: this
        REAL(MK), DIMENSION(:,:), POINTER :: d_theta
        INTEGER, INTENT(OUT)              :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_theta))THEN
           DEALLOCATE(d_theta)
        END IF
        
        ALLOCATE(d_theta(3,this%num_colloid))
        
        d_theta(:,:) = this%theta(:,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_theta


      SUBROUTINE colloid_get_omega(this,d_omega,stat_info)
        !---------------------------------------
        ! Return rotating velocity of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_omega
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_omega)) THEN           
           DEALLOCATE(d_omega)
        END IF
        
        ALLOCATE(d_omega(3,this%num_colloid))
        
        d_omega(1:3,1:this%num_colloid) = &
             this%omega(1:3,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_omega
      
      
      SUBROUTINE colloid_get_torque(this,d_torque,stat_info)
        !---------------------------------------
        ! Return torque of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_torque
        INTEGER, INTENT(OUT)                    :: stat_info
        
        
        stat_info = 0
        
        IF(ASSOCIATED(d_torque)) THEN           
           DEALLOCATE(d_torque)
        END IF
        
        ALLOCATE(d_torque(3,this%num_colloid))
        
        d_torque(1:3,1:this%num_colloid) = &
             this%torque(1:3,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_torque
      
      
      SUBROUTINE colloid_get_num_physical_part(this,&
           d_num_physical_part, stat_info)
        !---------------------------------------
        !  Return the number of physical particles
        !  which constitute colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, DIMENSION(:), POINTER  :: d_num_physical_part        
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_num_physical_part))THEN
           DEALLOCATE(d_num_physical_part)
        END IF

        ALLOCATE(d_num_physical_part(this%num_colloid))
       
        d_num_physical_part(:) = &
             this%num_physical_part(1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_num_physical_part
      
      
      INTEGER FUNCTION colloid_get_num_physical_part_tot(this,stat_info)
        !---------------------------------------
        !  Return total number of physical particles
        !  which constitute colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
      
        stat_info = 0 
        
        colloid_get_num_physical_part_tot = &
             this%num_physical_part_tot        
        
        RETURN       
        
      END FUNCTION colloid_get_num_physical_part_tot
      
      
      SUBROUTINE colloid_get_num_numerical_part(this,&
           d_num_numerical_part, stat_info)
        !---------------------------------------
        ! Return the number of numerical particles
        ! which constitute colloids.
        ! They are the outer layer of colloids.
        ! The width of layer is D*cut_off.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, DIMENSION(:), POINTER  :: d_num_numerical_part
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_num_numerical_part))THEN
           DEALLOCATE(d_num_numerical_part)
        END IF

        ALLOCATE(d_num_numerical_part(this%num_colloid))
      
        d_num_numerical_part(:) = &
             this%num_numerical_part(1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_num_numerical_part
      
      
      INTEGER FUNCTION colloid_get_num_numerical_part_tot(this,stat_info)
        !---------------------------------------
        !  Return total number of numerical particles
        !  which constitute colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
      
        stat_info = 0
        
        colloid_get_num_numerical_part_tot = &
             this%num_numerical_part_tot        
        
        RETURN       
        
      END FUNCTION colloid_get_num_numerical_part_tot
      

      SUBROUTINE colloid_get_f(this,d_f,stat_info)
        !--------------------------------------------
        ! Return translating accelaration of colloids.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_f
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_f)) THEN
           DEALLOCATE(d_f)
        END IF

        ALLOCATE(d_f(this%num_dim,this%num_colloid))
        
        d_f(:,:) = &
             this%f(1:this%num_dim,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_f
      

      REAL(MK) FUNCTION colloid_get_fa_min(this,stat_info)
        
        TYPE(Colloid), INTENT(IN)               :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        colloid_get_fa_min = this%fa_min
        
        RETURN       
        
      END FUNCTION  colloid_get_fa_min


      REAL(MK) FUNCTION colloid_get_fa_max(this,stat_info)
        
        TYPE(Colloid), INTENT(IN)               :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        colloid_get_fa_max = this%fa_max
        
        RETURN       
        
      END FUNCTION  colloid_get_fa_max
      
      
      REAL(MK) FUNCTION colloid_get_dt_f(this,stat_info)
        
        TYPE(Colloid), INTENT(IN)               :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        colloid_get_dt_f = this%dt_f
        
        RETURN       
        
      END FUNCTION  colloid_get_dt_f
      
      
      SUBROUTINE colloid_get_alpha(this,d_alpha,stat_info)
        !---------------------------------------
        ! Return rotating accleration colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_alpha
        INTEGER, INTENT(OUT)                    :: stat_info
        
        
        stat_info = 0
        
        IF(ASSOCIATED(d_alpha)) THEN           
           DEALLOCATE(d_alpha)
        END IF
        
        ALLOCATE(d_alpha(3,this%num_colloid))
        
        d_alpha(1:3,1:this%num_colloid) = &
             this%alpha(1:3,1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_alpha
      

      SUBROUTINE colloid_get_k_energy(this,d_k,stat_info)
        !--------------------------------------------
        ! Return the kinetic energy of colloids.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), POINTER          :: d_k
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_k)) THEN
           DEALLOCATE(d_k)
        END IF

        ALLOCATE(d_k(this%num_colloid))
        
        d_k(:) = &
             this%k_energy(1:this%num_colloid)
        
        RETURN       
        
      END SUBROUTINE colloid_get_k_energy
      
      
      SUBROUTINE colloid_get_mom(this,d_mom,stat_info)
        !--------------------------------------------
        ! Return the momentum  of colloids.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_mom
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_mom)) THEN           
           DEALLOCATE(d_mom)
        END IF

        ALLOCATE(d_mom(this%num_dim,this%num_colloid))
        
        d_mom(:,:) = &
             this%mom(1:this%num_dim,1:this%num_colloid)
        
        RETURN
        
      END SUBROUTINE colloid_get_mom
      
     
      REAL(MK) FUNCTION  colloid_get_k_energy_tot(this,stat_info)
        !--------------------------------------------
        ! Return the kinetic energy of all colloids.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        colloid_get_k_energy_tot = &
             this%k_energy_tot

        RETURN
        
      END FUNCTION colloid_get_k_energy_tot

      
      SUBROUTINE colloid_get_mom_tot(this,d_mom_tot,stat_info)
        !--------------------------------------------
        ! Return the momentum  of colloids.
        !--------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), POINTER         :: d_mom_tot
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_mom_tot)) THEN           
           DEALLOCATE(d_mom_tot)
        END IF

        ALLOCATE(d_mom_tot(this%num_dim))
        
        d_mom_tot(1:this%num_dim) = &
             this%mom_tot(1:this%num_dim)
        
        RETURN       
        
      END SUBROUTINE colloid_get_mom_tot

      
      REAL(MK) FUNCTION colloid_get_dout(this, stat_info)
        
        !---------------------------------------
        ! Get the minimal distance of a fluid
        ! particle from the surface of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        colloid_get_dout = this%dout
        
        RETURN
        
      END FUNCTION colloid_get_dout
      
      
      REAL(MK) FUNCTION colloid_get_din(this, stat_info)
        
        !---------------------------------------
        ! Get the maximum distance of a boundary
        ! particle from the surface of colloids.
        !---------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        stat_info = 0
        
        colloid_get_din = this%din
        
        RETURN
        
      END FUNCTION colloid_get_din
      
      
      REAL(MK) FUNCTION colloid_get_eta(this, stat_info)
        !---------------------------------------
        !   Return dynamic viscosity.
        !---------------------------------------

        TYPE(Colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        colloid_get_eta = this%eta
        
        RETURN
        
      END FUNCTION colloid_get_eta
      
      
      SUBROUTINE colloid_get_particle_v(this, &
           dim,b_x,b_v,b_sid,stat_info)
        
        !----------------------------------------------------
        ! Return the velocity of colloids.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        INTEGER, INTENT(IN)                     :: dim
        REAL(MK), DIMENSION(:), INTENT(IN)      :: b_x
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: b_v
        INTEGER, INTENT(IN)                     :: b_sid
        INTEGER, INTENT(OUT)                    :: stat_info
        
	!----------------------------------------------------
        ! Local variables start here :
	!----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        REAL(MK)                                :: d_max
        REAL(MK), DIMENSION(3)                  :: r_px
        REAL(MK), DIMENSION(3)                  :: r_pv
        REAL(MK), DIMENSION(3)                  :: x_coll
        REAL(MK), DIMENSION(3)                  :: v_coll
        
        REAL(MK), DIMENSION(:), POINTER         :: length        
        REAL(MK), DIMENSION(:,:), POINTER       :: shear_v
        REAL(MK), DIMENSION(:,:), POINTER       :: shear_length
        INTEGER                                 :: i, k
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(length)
        NULLIFY(shear_v)
        NULLIFY(shear_length)
        
        IF ( dim /= this%num_dim ) THEN
           PRINT *, "colloid_get_partcle_v : ", &
                "dim /= num_dim !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Assign zero velocity initially.
        !----------------------------------------------------
        
        b_v(1:dim) = 0.0_MK
        
        ALLOCATE(length(1:dim))
        
        length(:) = &
             this%max_phys(:) - this%min_phys(:)
        
        
        CALL boundary_get_shear_v(this%boundary, &
             shear_v,stat_info_sub)
        CALL boundary_get_shear_length(this%boundary, &
             shear_length,stat_info_sub)
      
        !----------------------------------------------------
        ! Check the shape of this colloid and decide
        ! d_max, which is used for checking if the boundary
        ! particle is too far away from the colloid's center,
        ! i.e., colloid crossed boundary alread.
        !----------------------------------------------------
        
        SELECT CASE ( this%shape(b_sid) ) 
           
        CASE (1) 
           
           d_max = this%radius(1,b_sid) + this%dout
           
        CASE (2) 
           
           d_max = this%radius(1,b_sid) 
           
           IF ( this%radius(2,b_sid) > this%radius(1,b_sid) ) THEN
              
              d_max = this%radius(2,b_sid)
              
           END IF
           
           d_max = d_max + this%dout
           
        CASE (3) 
           
           d_max = this%radius(1,b_sid) + &
                this%radius(2,b_sid) + this%dout
           
        END SELECT
        
        
        !----------------------------------------------------
        ! If colloid sid is translating, we assign
        ! velocity of translation to boundary particle.
        !----------------------------------------------------
        
        IF ( this%translate ) THEN
           
           v_coll(1:dim) = this%v(1:dim,b_sid)
           x_coll(1:dim) = this%x(1:dim,b_sid)
           
           r_px(1:dim) = &
                b_x(1:dim) - x_coll(1:dim)
           
           !-------------------------------------------------
           ! Loop over all dimensions.
           !-------------------------------------------------
           
           DO i = 1, dim
              
              !----------------------------------------------
              ! If the boundary particle is too far away
              ! from center of the colloid,
              ! it left the boundary from the max side
              ! and entered in min side already. 
              ! Then we have to use the image
              ! of the center.
              !----------------------------------------------
              
              IF( r_px(i) > d_max ) THEN
                 
                 !-------------------------------------------
                 ! Check for Lees-Edwards boundary.
                 !-------------------------------------------
                 
                 SELECT CASE (this%bcdef (2*i-1) )
                    
                 CASE ( ppm_param_bcdef_LE )
                    
                    DO k = 1, dim
                       
                       IF ( k /= i ) THEN
                          
                          v_coll(k) = v_coll(k) + &
                               shear_v(k,2*i) - shear_v(k,2*i-1)
                          
                       END IF
                       
                    END DO ! k = 1, dim
                    
                 END SELECT ! bcdef(2*i-1)
                 
                 !-------------------------------------------
                 ! If the boundary particle is too far away
                 ! from center of  the colloid,
                 ! it left the boundary from the min side
                 ! and entered from max side already.
                 ! Then we have to use the image of the center.
                 !-------------------------------------------
                 
              ELSE IF( -r_px(i) > d_max ) THEN
                 
                 SELECT CASE ( this%bcdef(2*i) )
                    
                 CASE ( ppm_param_bcdef_LE )
                    
                    DO k = 1, dim
                       
                       IF ( k /= i ) THEN
                          
                          v_coll(k) = v_coll(k) + &
                               shear_v(k,2*i-1) - shear_v(k,2*i)
                          
                       END IF
                       
                    END DO ! k = 1, dim
                    
                 END SELECT ! bcdef(2*i)
                 
              END IF ! r_px
              
           END DO ! i = 1, dim
           
           !-------------------------------------------------
           ! Add the translation velocity of nearst
           ! image for the center of the colloid.
           !-------------------------------------------------
           
           b_v(1:dim) =  &
                b_v(1:dim) + v_coll(1:dim)
           
        END IF
        
        !----------------------------------------------------
        ! If colloid sid is rotating, we assign velocity
        ! of rotation.
        !----------------------------------------------------
        
        IF ( this%rotate ) THEN
           
           x_coll(1:dim) = this%x(1:dim,b_sid)
           
           r_px(1:3)  = 0.0_MK              
           r_px(1:dim) = &
                b_x(1:dim) - x_coll(1:dim)
           
           !-------------------------------------------------
           ! Loop over all dimensions.
           !-------------------------------------------------
           
           DO i = 1, dim
              
              !----------------------------------------------
              ! If the boundary particle is too far away
              ! from center of the colloid, it left the
              ! boundary from the max side already.
              ! Then we have to use the image of the center.
              !----------------------------------------------
              
              IF( r_px(i) > d_max ) THEN
                 
                 SELECT CASE ( this%bcdef (2*i-1) )
                    
                 CASE  ( ppm_param_bcdef_periodic )  
                    
                    x_coll(i) = x_coll(i) + length(i)
                    
                 CASE ( ppm_param_bcdef_LE )
                    
                    DO k = 1, dim
                       
                       IF ( k== i ) THEN
                          
                          x_coll(k) = x_coll(k) + length(k)
                          
                       ELSE
                          
                          x_coll(k) = &
                               MODULO(x_coll(k) - shear_length(k,2*i-1), &
                               length(k) )
                          
                       END IF
                       
                    END DO ! k = 1, dim
                    
                 END SELECT ! bcdef(2*i-1)
                 
                 !-------------------------------------------
                 ! If the boundary particle is too far away
                 ! from center of the colloid, it crossed 
                 ! the boundary from the max side already.
                 ! Then we have to use the image of the center.
                 !-------------------------------------------
                 
              ELSE IF( -r_px(i) > d_max ) THEN
                 
                 SELECT CASE ( this%bcdef(2*i) )
                    
                 CASE ( ppm_param_bcdef_periodic )
                    
                    x_coll(i) = x_coll(i) - length(i)
                    
                 CASE ( ppm_param_bcdef_LE )
                    
                    DO k = 1, dim
                       
                       IF ( k== i ) THEN
                          
                          x_coll(k) = x_coll(k) - length(k)
                          
                       ELSE
                          
                          x_coll(k) = &
                               MODULO(x_coll(k) - shear_length(k,2*i), &
                               length(k) )
                          
                       END IF
                       
                    END DO ! k = 1, dim
                    
                 END SELECT ! bcdef(2*i)
                 
              END IF ! r_px(i)
              
           END DO ! i = 1, dim
           
           !----------------------------------------------
           ! Calculate the displacement of the boundary 
           ! particle from its center, then use angular
           ! velocity  of the colloid to get boundary 
           ! particles' equavilent velocity of translation.
           !----------------------------------------------
           
           r_px(1:dim) = b_x(1:dim) - x_coll(1:dim)
           
           CALL tool_cross_product(this%tool,&
                this%omega(1:3,b_sid),r_px(1:3),&
                r_pv(1:3),stat_info_sub)
           
           b_v(1:dim) = b_v(1:dim) + &
                r_pv(1:dim)
           
        END IF ! rotate(b_sid)
        
9999    CONTINUE
        
        IF(ASSOCIATED(length)) THEN
           DEALLOCATE(length)
        END IF
        
        IF(ASSOCIATED(shear_v)) THEN
           DEALLOCATE(shear_v)
        END IF
        
        IF(ASSOCIATED(shear_length)) THEN
           DEALLOCATE(shear_length)
        END IF
        
        RETURN       
        
      END SUBROUTINE colloid_get_particle_v
      
