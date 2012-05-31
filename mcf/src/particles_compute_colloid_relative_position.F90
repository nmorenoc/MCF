      SUBROUTINE particles_compute_colloid_relative_position(this,&
           stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_compute_colloid_relative_position
        !----------------------------------------------------
        !
        ! Purpose     : Compute the colloid boundary particle's
        !               relative position to the colloid center
        !               after rotation.
        !               
        !
        ! Reference   :
        !
        ! Remark      : Colloid are modelled as rigid body,
        !               therefore, relative position of a 
        !               colloid boundary particle to its
        !               center is calculated using rotation
        !               matrix.
        !
        ! Revision    : V0.2  17.11.2011, reformulated.
        !
        !               V0.1  23.08.2010, original version.
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
        !
        ! this           : an object of Particles Class.
        ! stat_info      : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------

        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim, i
        INTEGER                                 :: ip, sid
        TYPE(Colloid),POINTER                   :: colloids
        LOGICAL                                 :: translate,rotate
        REAL(MK), DIMENSION(3)                  :: coll_x, coll_v
        REAL(MK), DIMENSION(:,:,:), POINTER     :: coll_rot_matrix
        REAL(MK), DIMENSION(3)                  :: tx, rx_a, rx_b
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        rotate        = .FALSE.
        
        dim = physics_get_num_dim(this%phys,stat_info_sub)
        
        NULLIFY(colloids)
        CALL physics_get_colloid(this%phys,colloids,stat_info_sub)
        translate    = &
             colloid_get_translate(colloids,stat_info_sub)

        rotate    = &
             colloid_get_rotate(colloids,stat_info_sub)
        NULLIFY(coll_rot_matrix)
        
#ifdef __PARTICLES_POSITION_FIXED
#else
        
        IF ( translate .OR. rotate ) THEN
           
           !-------------------------------------------------
           ! Get rotation matrix of colloids.
           !-------------------------------------------------
           
           IF ( rotate ) THEN
              
              CALL colloid_get_rotation_matrix(colloids,&
                   coll_rot_matrix,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "particles_integrate_colloid_position : " ,&
                      "Getting colloid rot_matrix failed!"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF
           
           !-------------------------------------------------
           ! Loop over all colloid boundary particles.
           !-------------------------------------------------
           
           DO i = 1, this%num_part_colloid
              
              !----------------------------------------------
              ! Get index of this boundary particle
              ! and its species ID.
              !----------------------------------------------
              
              ip  = this%part_colloid_list(1,i)
              sid = this%part_colloid_list(2,i)
              
              !----------------------------------------------
              ! Get the nearest image of colloid center.
              ! Note that coll_x must be from previous
              ! step, is not updated yet!
              !----------------------------------------------
              
              CALL colloid_in_nearest_image(colloids, &
                   this%x(1:dim,ip),sid,coll_x,rx_a,&
                   coll_v,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "particles_compute_colloid_relative_position: " ,&
                      "Getting colloid nearest image failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! For rotation:
              ! get new relative position according to 
              ! rotation motion.
              !
              ! For non-rotation:
              ! keep the old relative position.
              !----------------------------------------------
              IF ( rotate ) THEN
              
                 rx_b(1:dim) = &
                      MATMUL(coll_rot_matrix(1:dim,1:dim,sid),rx_a(1:dim))
                 this%x(1:dim,ip) = rx_b(1:dim) 

              ELSE

                 this%x(1:dim,ip) = rx_a(1:dim)
                 
              END IF
              
           END DO ! i =1, num_part_colloid
           
        END IF ! translate .OR. rotate
           
#endif        
        
9999    CONTINUE

        IF ( ASSOCIATED(coll_rot_matrix) ) THEN
           DEALLOCATE(coll_rot_matrix)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_compute_colloid_relative_position
      
      
