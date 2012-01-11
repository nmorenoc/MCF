      SUBROUTINE particles_adjust_particles(this,num,stat_info)
	!----------------------------------------------------
     	! Subroutine  : particles_adjust_particles
	!----------------------------------------------------
      	!
      	! Purpose     : Adjust particles positions and 
        !               velocities according to boundary
        !               conditions.
        !                
      	!	 	      	 
        ! Reference   :
        !              
      	! Remark      : In case of particle crossing
        !               boundary line,
        !
        !               1  for periodic boundary, reinsert
        !                  particle on the other side with
        !                  same velocity.
        !
        !               3  for symmetry boundary, reflect
        !                  fluid particle back with 
        !                  opposite normal velocity;
        !                  freeze colloid boundary particle.
        !
        !               7  for wall using symmetry/mirror or
        !               9  for wall using solid particles,
        !                  freez fluid/colloid particle
        !                  on the border with zero velocity
        !                  component normal to the wall.
        !
        !               10 for Lees-Edwards boundary,
        !                  consider the shear relative
        !                  velocity and displacement
        !                  between layers of boxes.
        !
        !               According to ppm_impose_part_bc.F90 routine
        !               IF statements matters and should be seperated
        !               (no ELSEIFs), all due to round off errors.
      	!
        ! Revisions   : V0.3 18.11 2009, including
        !               Lees-Edwards boundaries.
        !
        !               V0.2 08.07 2009,
        !               check again the work flow is correct,
        !               supply with more comments.
        !                 
        !               V0.1 11.02 2009, original version.
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
        ! num        : number of particles needed to be updated,
        !              i.e. first 'num' particles in this%x are 
        !              operated.
        ! stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)  :: this
        INTEGER, INTENT(IN)             :: num
        INTEGER, INTENT(OUT)	        :: stat_info
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(:), POINTER :: min_phys
        REAL(MK), DIMENSION(:), POINTER :: max_phys
        REAL(MK), DIMENSION(:), POINTER :: length
        REAL(MK), DIMENSION(:), POINTER :: min_phys_t
        REAL(MK), DIMENSION(:), POINTER :: max_phys_t      
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        TYPE(Boundary), POINTER         :: tboundary
        REAL(MK),DIMENSION(:,:),POINTER :: shear_length
        REAL(MK),DIMENSION(:,:),POINTER :: shear_v
        INTEGER	 		        :: i,j, k
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(min_phys)
        NULLIFY(max_phys)
        NULLIFY(min_phys_t)
        NULLIFY(max_phys_t)   
        NULLIFY(length)
        NULLIFY(bcdef)
        NULLIFY(tboundary)
        NULLIFY(shear_length)
        NULLIFY(shear_v)
                
        !----------------------------------------------------
        ! Check if num is in the range of particles' number,
        ! we are supposed to check only the real particles.
        !----------------------------------------------------
        
        IF( num > this%num_part_real ) THEN
           PRINT *, "particles_adjust_particles :", &
                "Num > number of real particles !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Get physics quantities including boundary.
        !----------------------------------------------------
        
        num_dim = physics_get_num_dim(this%phys,stat_info_sub)
        CALL physics_get_min_phys(this%phys,min_phys,stat_info_sub)
        CALL physics_get_max_phys(this%phys,max_phys,stat_info_sub)
        CALL physics_get_length(this%phys,length,stat_info_sub)
        CALL physics_get_min_phys_t(this%phys,min_phys_t,stat_info_sub)
        CALL physics_get_max_phys_t(this%phys,max_phys_t,stat_info_sub)
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        CALL boundary_get_shear_length(tboundary,shear_length,stat_info_sub)
        CALL boundary_get_shear_v(tboundary,shear_v,stat_info_sub)
        
        !----------------------------------------------------
        ! Loop over the first num particles.
        !----------------------------------------------------
        
        DO j = 1, num
	   
           !-------------------------------------------------
           ! Check each dimension.
           !-------------------------------------------------
           
           DO i = 1, num_dim
              
              !----------------------------------------------
              ! Check if particle j is out of lower boundary.
              ! with '<' instead of '<='.
              !----------------------------------------------
              
              IF ( this%x(i,j) <  min_phys(i) ) THEN
                 
                 SELECT CASE( bcdef(2*i-1) )
                    
                 CASE ( ppm_param_bcdef_periodic ) 
                    
                    !----------------------------------------
                    ! Apply to each kind of particle,
                    ! reinsert it back from the other side.
                    !----------------------------------------
                    
                    this%x(i,j) = this%x(i,j) + length(i)
                    
                 CASE ( ppm_param_bcdef_symmetry ) 
                    
                    IF ( this%id(this%sid_idx,j) == 0 )THEN
                       
                       !-------------------------------------
                       ! Apply to fluid particle, reflect it
                       ! back, give opposite normal velocity.
                       !-------------------------------------

                       this%x(i,j) = &
                            2.0_MK * min_phys(i) - this%x(i,j)
                       
                       this%v(i,j) = - this%v(i,j)
                       
                    ELSE IF ( this%id(this%sid_idx,j) > 0 )THEN
                       
                       !-------------------------------------
                       ! freeze colloid particle locally.
                       !-------------------------------------
                       
                       this%x(i,j) = min_phys(i)
                       
                       this%v(i,j) = 0.0_MK
                       
                    END IF
                    
                 CASE ( ppm_param_bcdef_wall_sym )
                    
                    !----------------------------------------
                    ! Apply to all particles, freeze its
                    ! normal direction.
                    ! Not ghost particle
                    !----------------------------------------
                    
                    !this%x(i,j) = min_phys(i)
                    
                    !this%v(i,j) = 0.0_MK
                    
                 CASE ( ppm_param_bcdef_wall_solid )
                    
                    IF ( this%id(this%sid_idx,j) >= 0 ) THEN
                       
                       !-------------------------------------
                       ! Apply to fluid and colloid boundary
                       ! particle, freeze its normal direction.
                       !-------------------------------------
                       
                       !this%x(i,j) = min_phys(i)
                       
                       !this%v(i,j) = 0.0_MK
                       
                    ELSE IF ( this%id(this%sid_idx,j) < 0 .AND. &
                         this%x(i,j) < min_phys_t(i) ) THEN
                       
                       !-------------------------------------
                       ! Apply to a solid wall particle,
                       ! freeze its normal direction.
                       !-------------------------------------
             
                       this%x(i,j) = min_phys_t(i)
                       
                       this%v(i,j) = 0.0_MK                      
                       
                    END IF
                    
                 CASE ( ppm_param_bcdef_LE )
                    
                    IF ( this%id(this%sid_idx,j) >= 0 ) THEN
                    
                       !-------------------------------------
                       ! Apply to fluid/colloid particle.
                       !-------------------------------------
                   
                       DO k = 1, num_dim
                          
                          IF ( k == i ) THEN
                             
                             this%x(k,j) = this%x(k,j) + length(k)
                             
                          ELSE
                             
                             this%x(k,j) = &
                                  MODULO( this%x(k,j) - shear_length(k,2*i-1), &
                                  length(k) )
                             
                             this%v(k,j) = this%v(k,j) + &
                                  ( shear_v(k,2*i) - shear_v(k,2*i-1) )
                             
                          END IF
                          
                       END DO ! k = 1, num_dim
                       
                    END IF ! id(sid_idx,j) >=0
                    
                    EXIT ! Found one dimension, quit.

                 END SELECT ! bcdef(2*i-1)                 
                 
              END IF ! x < min_phys
              
              !-------------------------------------------
              ! Check if j is out of upper boundary.
              ! With '>=' instead of '>'.
              !-------------------------------------------
              
              IF ( this%x(i,j) >=  max_phys(i) ) THEN
                 
                 SELECT CASE ( bcdef(2*i) )
                    
                 CASE ( ppm_param_bcdef_periodic )
                    
                    !----------------------------------------
                    ! Apply to each kind of particle,
                    ! reinsert it from other side.
                    !----------------------------------------
                    
                    this%x(i,j) = this%x(i,j) - length(i)
                    
                 CASE ( ppm_param_bcdef_symmetry ) 
                    
                    IF ( this%id(this%sid_idx,j) == 0 )THEN
                       
                       !-------------------------------------
                       ! Apply to fluid particle, reflect it
                       ! back, give opposite normal velocity.
                       !-------------------------------------
                       
                       this%x(i,j) = &
                            2.0_MK * max_phys(i) - this%x(i,j)
                       
                       this%v(i,j) = - this%v(i,j)
                       
                    ELSE IF ( this%id(this%sid_idx,j) > 0 )THEN
                       
                       !-------------------------------------
                       ! Apply to colloid particle, freeze 
                       ! its normal direction.
                       !-------------------------------------
                       
                       this%x(i,j) = &
                            max_phys(i) - mcf_machine_zero
                       
                       this%v(i,j) = 0.0_MK
                       
                    END IF
                    
                 CASE ( ppm_param_bcdef_wall_sym )
                    
                    !----------------------------------------
                    ! Apply to all particles, freeze its
                    ! normal direction.
                    ! Not ghost particle
                    !----------------------------------------
                    
                    !this%x(i,j) = &
                    !     max_phys(i) - mcf_machine_zero
                    
                    !this%v(i,j) = 0.0_MK
                    
                    
                 CASE ( ppm_param_bcdef_wall_solid )
                    
                    IF ( this%id(this%sid_idx,j) >= 0 ) THEN
                       
                       !-------------------------------------
                       ! Apply to fluid and colloid boundary
                       ! particle, freeze its normal direction.
                       !-------------------------------------
                       
                     !  this%x(i,j) = &
                     !       max_phys(i) - mcf_machine_zero
                       
                     !  this%v(i,j) = 0.0_MK
                       
                    ELSE IF ( this%id(this%sid_idx,j) < 0 .AND. &
                         this%x(i,j) > max_phys_t(i) ) THEN
                       
                       !-------------------------------------
                       ! Apply to a solid wall particle,
                       ! freeze its normal direction.
                       !-------------------------------------
                       
                       this%x(i,j) = &
                            max_phys_t(i) - mcf_machine_zero
                       
                       this%v(i,j) = 0.0_MK   
                       
                    END IF
                    
                 CASE ( ppm_param_bcdef_LE ) 
                    
                    !----------------------------------------
                    ! Apply to a fluid/colloid particle.
                    !----------------------------------------
                    
                    IF ( this%id(this%sid_idx,j) >= 0 ) THEN
                       
                       DO k = 1, num_dim
                          
                          IF ( k == i ) THEN
                             
                             this%x(k,j) = this%x(k,j) - length(k)
                             
                          ELSE
                             
                             this%x(k,j) = &
                                  MODULO( this%x(k,j) - shear_length(k,2*i), &
                                  length(k) )
                             
                             this%v(k,j) = this%v(k,j) + &
                                  ( shear_v(k,2*i-1) - shear_v(k,2*i) )
                             
                          END IF ! k==i
                          
                       END DO ! k = 1, num_dim
                       
                    END IF ! id(sid_idx,j) > = 0
                    
                    EXIT ! Found one dimension, quit.
                    
                 END SELECT ! bcdef(2*i)
                 
              END IF ! x > max_phys
              
           END DO ! i = 1, num_dim
           
        END DO ! j = 1, num of particles
        
        
9999    CONTINUE	
        
        IF(ASSOCIATED(min_phys)) THEN
           DEALLOCATE(min_phys) 
        END IF
        
        IF(ASSOCIATED(max_phys)) THEN
           DEALLOCATE(max_phys) 
        END IF
        
        IF(ASSOCIATED(length)) THEN
           DEALLOCATE(length) 
        END IF
        
        IF(ASSOCIATED(min_phys_t)) THEN
           DEALLOCATE(min_phys_t) 
        END IF
        
        IF(ASSOCIATED(max_phys_t)) THEN
           DEALLOCATE(max_phys_t) 
        END IF
        
        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef) 
        END IF
        
        IF(ASSOCIATED(shear_length)) THEN
           DEALLOCATE(shear_length) 
        END IF
        
        IF(ASSOCIATED(shear_v)) THEN
           DEALLOCATE(shear_v) 
        END IF
        
        
        RETURN
        
      END SUBROUTINE particles_adjust_particles
