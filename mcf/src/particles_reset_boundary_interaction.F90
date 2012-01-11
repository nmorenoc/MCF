      SUBROUTINE particles_reset_boundary_interaction(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_reset_boundary_interaction
        !----------------------------------------------------
        !
        ! Purpose     : Assign solid wall particles with
        !               zero force.
        !
        ! Refernece   :
        !
        ! Remark      : 
        !      
        !
        ! Revision    : V0.1 22.11 20101, original version.
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
        ! this       : an object of Particles Class.
        ! stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------
        
        INTEGER                                 :: dim
        INTEGER                                 :: i,ip
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        dim       = this%num_dim
        
        !----------------------------------------------------
        ! Assign particles which constitute solid walls 
        ! with zero force.
        !----------------------------------------------------
        
#if 0
        DO i =1, this%num_part_wall_solid_real
           
           ip  = this%part_wall_solid_real_list(1,i)
           
           this%f(1:dim,ip) = 0.0_MK
           
        END DO
#endif   
        
9999    CONTINUE
        
        
        RETURN
        
      END SUBROUTINE particles_reset_boundary_interaction
      
      
      
