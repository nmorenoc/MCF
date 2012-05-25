      SUBROUTINE particles_compute_colloid_absolute_position(this,&
           stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_compute_colloid_absolute_position
        !----------------------------------------------------
        !
        ! Purpose     : Compute the colloid boundary particle's
        !               absolute position after the colloid center
        !               translated.
        !               
        !
        ! Reference   :
        !
        ! Remark      : Colloid are modelled as rigid body.
        !
        ! Revision    : V0.1  17.11.2011, original version.
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
        LOGICAL                                 :: translate, rotate
        REAL(MK), POINTER, DIMENSION(:,:)       :: coll_x
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        rotate        = .FALSE.
        
        dim = physics_get_num_dim(this%phys,stat_info_sub)
        
        NULLIFY(colloids)
        NULLIFY(coll_x)
        
        CALL physics_get_colloid(this%phys,colloids,stat_info_sub)
        translate    = &
             colloid_get_translate(colloids,stat_info_sub)
        rotate       = &
             colloid_get_rotate(colloids,stat_info_sub)
        
        CALL colloid_get_x(colloids,coll_x,stat_info_sub)
        
#ifdef __PARTICLES_POSITION_FIXED
#else
        
        IF ( translate .OR. rotate ) THEN
           
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
              
              this%x(1:dim,ip) = this%x(1:dim,ip) + coll_x(1:dim,sid)
              
           END DO ! i =1, num_part_colloid
           
        END IF ! translate OR rotate
        
#endif        
        
9999    CONTINUE

        IF ( ASSOCIATED(coll_x) ) THEN
           DEALLOCATE(coll_x)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_compute_colloid_absolute_position
      
