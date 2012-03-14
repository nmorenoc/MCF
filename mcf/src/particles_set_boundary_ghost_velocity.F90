      SUBROUTINE particles_set_boundary_ghost_velocity(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_set_boundary_ghost_velocity
        !----------------------------------------------------
        !
        ! Purpose     : Set boundary ghost particles' 
        !               velocity to have prescribed 
        !               boundary conditions.
        !                 
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     : 
        !                
        !
        ! Revisions   : V0.2 14.10 2010, change <= to <
        !               for min_phys side. 
        !
        !               V0.1 07.12 2009, original version.
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
        ! Physics parameters :
        !----------------------------------------------------

        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(:), POINTER :: min_phys
        REAL(MK), DIMENSION(:), POINTER :: max_phys

        !----------------------------------------------------
        ! Boundary parameters :
     	!----------------------------------------------------
        
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        TYPE(Boundary), POINTER         :: tboundary
        REAL(MK),DIMENSION(:,:),POINTER :: shear_v
        INTEGER                         :: num_sym
        INTEGER                         :: num_wall_sym
        INTEGER                         :: num_wall_solid
        INTEGER                         :: num_le
        INTEGER                         :: i, j, k
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        NULLIFY(min_phys)
        NULLIFY(max_phys)

        NULLIFY(bcdef)
        NULLIFY(tboundary)        
        NULLIFY(shear_v)

        !----------------------------------------------------
        ! Physics_parameters :
        !----------------------------------------------------
        
        num_dim = this%num_dim

        CALL physics_get_min_phys(this%phys,min_phys,stat_info_sub)
        CALL physics_get_max_phys(this%phys,max_phys,stat_info_sub)
        
        !----------------------------------------------------
        ! Boundary parameters :
        !----------------------------------------------------
        
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)        
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        CALL boundary_get_shear_v(tboundary,shear_v,stat_info_sub)
        
        num_sym      = &
             boundary_get_num_sym(tboundary,stat_info_sub)
        num_wall_sym = &
             boundary_get_num_wall_sym(tboundary,stat_info_sub)
        num_wall_solid = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        num_le       = &
             boundary_get_num_le(tboundary,stat_info_sub)
        
        !----------------------------------------------------
        ! 1 For periodic boundaries:
        !
        ! we do nothing.
        !
        ! 3 For symmetry boundaries:
        !
        ! we assign mirror direction of velocity to ghost 
        ! particles, i.e., opposite velocity normal to boundary.
        !
        ! 7 For wall boundaries created by PPM using 
        !   symmetric/mirror particles:
        !
        ! we give opposite direction of velocity considering 
        ! the wall shear velocity, in order to have no slip 
        ! on the surface of walls;
        !
        ! 9 For wall using solid boundary particle:
        !
        ! it was  given correct wall velocity already before 
        !  creating ghost. we do nothing for it here.
        !
        ! 10 For Lees-Edwards boundaries:
        !
        ! we add relative shear velocity between 
        ! adjacent layers.
        !----------------------------------------------------
        
        
        !----------------------------------------------------
        ! If there is ghost particle making up the boundary
        ! condition, loop over each ghost particle.
        !----------------------------------------------------
        
        IF ( num_sym > 0 .OR. num_wall_sym > 0 .OR. &
             num_le > 0 ) THEN
           
           DO j = this%num_part_real + 1, this%num_part_all
              
              !----------------------------------------------
              ! Check each dimension at boundary and give 
              ! correct velocity to the particle 
              ! which is outside of boundary.
              !----------------------------------------------
              
              DO i = 1, num_dim
                 
                 !-------------------------------------------
                 ! For the one outside of lower boundary.
                 !-------------------------------------------
                 
                 IF( this%x(i,j) < min_phys(i) ) THEN
                    
                    SELECT CASE ( bcdef(2*i-1) )
                       
                    CASE ( ppm_param_bcdef_symmetry )
                       
                       !-------------------------------------
                       ! Reflect normal velocity.
                       !-------------------------------------
                       
                       this%v(i,j) = - this%v(i,j)
                       
                    CASE ( ppm_param_bcdef_wall_sym )
                       
                       !-------------------------------------
                       ! Give opposite velocity.
                       !-------------------------------------
                       
                       CALL boundary_noslip(tboundary,&
                            this%v(1:num_dim,j), 1-2*i,stat_info_sub)
                       
                    CASE ( ppm_param_bcdef_LE )
                       
                       !-------------------------------------
                       ! Plus the shear velocity.
                       !-------------------------------------
                       
                       DO k = 1, num_dim
                          
                          IF ( k /=i ) THEN
                             
                             this%v(k,j) = this%v(k,j) + &
                                  shear_v(k,2*i-1) - shear_v(k,2*i)
                             
                          END IF ! k/=i
                          
                       END DO ! k
                       
                    END SELECT ! bcdef(2*i-1)
                    
                    !----------------------------------------
                    ! For the one outside of upper boundary.
                    !----------------------------------------
                    
                 ELSE IF( this%x(i,j) >= max_phys(i) ) THEN
                    
                    SELECT CASE ( bcdef(2*i) )
                       
                    CASE ( ppm_param_bcdef_symmetry )
                       
                       !-------------------------------------
                       ! Reflect normal velocity.
                       !-------------------------------------
                       
                       this%v(i,j) = - this%v(i,j)
                       
                    CASE (ppm_param_bcdef_wall_sym )
                       
                       !-------------------------------------
                       ! Give opposite velocity.
                       !-------------------------------------
                       
                       CALL boundary_noslip(tboundary,&
                            this%v(1:num_dim,j), -2*i,stat_info_sub)
                       
                    CASE ( ppm_param_bcdef_LE )
                       
                       !-------------------------------------
                       ! Plus the shear velocity.
                       !-------------------------------------
                       
                       DO k = 1, num_dim
                          
                          IF ( k /=i ) THEN
                             
                             this%v(k,j) = this%v(k,j) + &
                                  shear_v(k,2*i) - shear_v(k,2*i-1)
                             
                          END IF ! k/=i
                          
                       END DO ! k
                       
                    END SELECT ! bcdef(2*i)
                    
                 END IF ! x(i,j) > = max_phys(i)
                 
              END DO ! i = 1, num_dim
              
           END DO ! j = real+1, all
           
        END IF

        
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
        
      END SUBROUTINE particles_set_boundary_ghost_velocity
      
      
