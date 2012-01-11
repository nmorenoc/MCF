      SUBROUTINE particles_reset_colloid_velocity(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_reset_colloid_velocity
        !----------------------------------------------------
        !
        ! Purpose     : Reset colloidal boundary particles
        !               velocity to zero.
        !
        ! Refernece   :
        !
        ! Remark      :
        !
        ! Revision    : V0.1 14.10.2009, original version.
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
        ! this       : an object of Particles Class.
        ! stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        
	!----------------------------------------------------
        ! Local variables start here :
	!----------------------------------------------------
        
      
        INTEGER                                 :: dim
        INTEGER                                 :: i,ip
        REAL(MK), DIMENSION(3)                  :: v_p
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info  = 0
        dim        = this%num_dim
        
        
        !----------------------------------------------------
        ! Loop over all colloidal boundary particles and
        ! assign them zero.
        !----------------------------------------------------
           
        DO i = 1, this%num_part_colloid
           
           !----------------------------------------------
           ! Get index of this boundary particle
           ! and its species ID.
           !----------------------------------------------
              
           ip  = this%part_colloid_list(1,i)
           
           this%v(1:dim,ip) = 0.0_MK
           
        END DO ! i = 1, num_part_colloid
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE particles_reset_colloid_velocity
      
      
