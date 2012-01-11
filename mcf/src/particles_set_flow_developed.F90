      SUBROUTINE particles_set_flow_developed(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_set_flow_developed
        !----------------------------------------------------
        !
        ! Purpose     : Set solvent particles velocity to
        !               the develop velocity. 
        !
        ! Remarks     : Currently, shear flow is available.
        !                
        !
        ! Revisions   :V0.1 3.12 2010, original version.
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
        INTEGER                         :: stat_info_sub
        INTEGER                         :: dim, i,j
        REAL(MK),DIMENSION(:),POINTER   :: min_phys
        REAL(MK),DIMENSION(:),POINTER   :: max_phys
        REAL(MK),DIMENSION(3)           :: length
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER,DIMENSION(:),POINTER    :: bcdef
        INTEGER                         :: num_wall_solid
        REAL(MK),DIMENSION(:,:),POINTER :: shear_v        

        !----------------------------------------------------
        ! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        NULLIFY(min_phys)
        NULLIFY(max_phys)
        NULLIFY(tboundary)
        NULLIFY(bcdef)
        NULLIFY(shear_v)

        !----------------------------------------------------
        ! Physics_parameters :
        !----------------------------------------------------
        
        dim = this%num_dim
        CALL physics_get_min_phys(this%phys,min_phys,stat_info_sub)
        CALL physics_get_max_phys(this%phys,max_phys,stat_info_sub)
        length(1:dim) = max_phys(1:dim) - min_phys(1:dim)
        
        !----------------------------------------------------
        ! Boundary parameters :
        !----------------------------------------------------
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        num_wall_solid = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        
        IF ( num_wall_solid > 0 ) THEN
           
           CALL boundary_get_shear_v(tboundary,shear_v,stat_info_sub)

           DO i = 1, dim
              
              IF ( bcdef(2*i) == ppm_param_bcdef_wall_solid ) THEN
                 
                 DO j = 1, this%num_part_real
                    
                    IF ( this%id(this%sid_idx,j) == &
                         mcf_particle_type_fluid ) THEN
                       
                       this%v(1:dim,j) = &
                            (shear_v(1:dim,2*i) - shear_v(1:dim,2*i-1)) * &
                            ( this%x(i,j) - min_phys(i)) / &
                            length(i) + shear_v(1:dim,2*i-1)
                       
                    END IF
                    
                 END DO ! j = 1 , num_part_real
                 
              END IF ! wall_solid
              
           END DO ! i = 1, dim
           
        END IF ! num_wall_solid > 0


        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Release all dynamics memories.
        !----------------------------------------------------
        
        IF(ASSOCIATED(min_phys)) THEN
           DEALLOCATE(min_phys)
        END IF
        
        IF(ASSOCIATED(max_phys)) THEN
           DEALLOCATE(max_phys)
        END IF
        
        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF

        IF(ASSOCIATED(shear_v)) THEN
           DEALLOCATE(shear_v)
        END IF

        RETURN
        
      END SUBROUTINE particles_set_flow_developed
      
