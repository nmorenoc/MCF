      SUBROUTINE particles_set_boundary_velocity(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_set_boundary_velocity
        !----------------------------------------------------
        !
        ! Purpose     : Set boundary condition boundary 
        !               particles velocity,
        !               to have correct boundary conditions.
        !
        ! Remarks     : 
        !                
        !
        ! Revisions   :V0.1 22.11 2010, original version.
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
        ! Argumetns :
      	!----------------------------------------------------
        
        TYPE(Particles),INTENT(INOUT)   :: this
        INTEGER,INTENT(OUT)             :: stat_info
        
        !----------------------------------------------------
        ! Physics, boundary parameters :
        !----------------------------------------------------

        INTEGER                         :: dim
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: num_wall_solid
        REAL(MK),DIMENSION(:,:),POINTER :: shear_v        
        INTEGER                         :: i,ip,sid
        
        INTEGER                         :: stat_info_sub

        !----------------------------------------------------
        ! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        NULLIFY(tboundary)        
        NULLIFY(shear_v)
        
        !----------------------------------------------------
        ! Physics_parameters :
        !----------------------------------------------------
        
        dim = this%num_dim
        
        !----------------------------------------------------
        ! Boundary parameters :
        !----------------------------------------------------
        
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        num_wall_solid = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        
        IF ( num_wall_solid > 0 ) THEN
           
           CALL boundary_get_shear_v(tboundary,shear_v,stat_info_sub)
        
           !-------------------------------------------------
           ! Assign particles which constitute solid walls 
           ! with solid walls' v.
           !-------------------------------------------------
           
           DO i =1, this%num_part_wall_solid_real
              
              ip  = this%part_wall_solid_real_list(1,i)
              sid = this%part_wall_solid_real_list(2,i)
              sid = ABS(sid)
              this%v(1:dim,ip) = shear_v(1:dim,sid)
              
           END DO
           
        END IF
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Release all dynamics memories.
        !----------------------------------------------------
        
        IF(ASSOCIATED(shear_v)) THEN
           DEALLOCATE(shear_v)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_set_boundary_velocity
      
      
