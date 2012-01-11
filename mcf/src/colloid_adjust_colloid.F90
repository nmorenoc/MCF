      SUBROUTINE colloid_adjust_colloid(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_adjust_colloid
        !----------------------------------------------------

        ! Purpose     : Adjust position, velocity of 
        !               centers of colloids as rigid bodies,
        !               accroding to boundary conditions.
        !
        ! Remark      : In case of particle crossing
        !               boundary line,
        !
        !               1  for periodic boundary, reinsert
        !                  center of particle on the other
        !                  side with same velocity.
        !
        !               3  for symmetry boundary, freeze
        !                  the center locally.
        !
        !               7  for wall using symmetry/mirror
        !                  freez it locally.
        !
        !               9  for wall using solid particles,
        !                  freez the center locally.
        !
        !               10 for Lees-Edwards boundary,
        !                  consider the shear relative 
        !                  velocity and displacement
        !                  between layer boxes.
        !          
        ! Note that this routine needs to be improved
        ! to stop colloid's surface away from wall not
        ! its center.
        !
        ! Revision    : V0.2 26.11.2009, including
        !               Lees-Edwards boundaries.
        !
        !               V0.1 15.07.2009, original version
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
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(OUT)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables start hear :
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(:), POINTER :: length
        REAL(MK),DIMENSION(:,:),POINTER :: shear_length
        REAL(MK),DIMENSION(:,:),POINTER :: shear_v
        INTEGER                         :: i,j,k
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(length)
        NULLIFY(shear_length)
        NULLIFY(shear_v)
        
        !----------------------------------------------------
        ! Get physics and boundary parametes.
        !----------------------------------------------------
        
        num_dim = this%num_dim
        
        ALLOCATE(length(num_dim))
        
        length(1:num_dim) = &
             this%max_phys(1:num_dim) - this%min_phys(1:num_dim)
     
        CALL boundary_get_shear_length(this%boundary,&
             shear_length,stat_info_sub)
        CALL boundary_get_shear_v(this%boundary,&
             shear_v,stat_info_sub)
        
        
        !----------------------------------------------------
        ! Loop over all colloids.
        !----------------------------------------------------
        
        DO j = 1, this%num_colloid
	   
           !-------------------------------------------------
           ! Check each dimension.
           !-------------------------------------------------
           
           DO i = 1, num_dim
              
              !----------------------------------------------
              ! Check if colloid j is out of lower boundary.
              !----------------------------------------------
              
              IF ( this%x(i,j) < this%min_phys(i) ) THEN
                 
                 SELECT CASE( this%bcdef(2*i-1) )
                    
                 CASE ( ppm_param_bcdef_periodic )
                    
                    !----------------------------------------
                    ! Reinsert it back from the other side.
                    !----------------------------------------
                    
                    this%x(i,j) = this%x(i,j) + &
                         length(i)
                    
                 CASE ( ppm_param_bcdef_symmetry )
                    
                    !----------------------------------------
                    ! Freeze its center locally.
                    !----------------------------------------
                    
                    this%x(i,j) = this%min_phys(i) 
                    
                    this%v(i,j) = 0.0_MK
                    
                 CASE ( ppm_param_bcdef_wall_sym ) 
                    
                    !----------------------------------------
                    ! Freeze its center locally.
                    !----------------------------------------
                    
                    this%x(i,j) = this%min_phys(i) 
                    
                    this%v(i,j) = 0.0_MK
                    
                 CASE ( ppm_param_bcdef_wall_solid )
                    
                    !----------------------------------------
                    ! Freeze its center locally.
                    !----------------------------------------
                    
                    this%x(i,j) = this%min_phys(i)
                    
                    this%v(i,j) = 0.0_MK
                    
                 CASE ( ppm_param_bcdef_LE ) 
                    
                    DO k = 1, num_dim
                       
                       IF ( k == i ) THEN
                          
                          this%x(k,j) = this%x(k,j) + length(k)
                          
                       ELSE
                          
                          this%x(k,j) = &
                               MODULO(this%x(k,j) - shear_length(k,2*i-1), &
                               length(k) )
                          
                          this%v(k,j) = this%v(k,j) + &
                               ( shear_v(k,2*i) - shear_v(k,2*i-1) )
                          
                       END IF ! k == i
                       
                    END DO ! k = 1, num_dim
                    
                    EXIT ! Found one dimension, quit.
   
                 END SELECT ! bcdef(2*i-1)
                 
              END IF ! x < min_phys
              
              !----------------------------------------------
              ! Check if j is out of upper boundary.
              !----------------------------------------------
              
              IF (this%x(i,j) >= this%max_phys(i) ) THEN
                 
                 SELECT CASE ( this%bcdef(2*i) ) 
                    
                 CASE ( ppm_param_bcdef_periodic )
                    
                    !----------------------------------------
                    ! Reinsert it back from the other side.
                    !----------------------------------------
                    
                    this%x(i,j) = this%x(i,j) - &
                         length(i)
                    
                 CASE ( ppm_param_bcdef_symmetry )
                    
                    !----------------------------------------
                    ! Freeze its center locally.
                    !----------------------------------------
                    
                    this%x(i,j) = &
                         this%max_phys(i) - mcf_machine_zero
                    
                    this%v(i,j) = 0.0_MK
                    
                 CASE ( ppm_param_bcdef_wall_sym ) 
                    
                    !----------------------------------------
                    ! Freeze its center locally.
                    !----------------------------------------
                    
                    this%x(i,j) = &
                         this%max_phys(i)  - mcf_machine_zero
                    
                    this%v(i,j) = 0.0_MK
                    
                 CASE ( ppm_param_bcdef_wall_solid )
                    
                    !----------------------------------------
                    ! Freeze its center locally.
                    !----------------------------------------
                    
                    this%x(i,j) = &
                         this%max_phys(i) - mcf_machine_zero
                    
                    this%v(i,j) = 0.0_MK
                    
                 CASE ( ppm_param_bcdef_LE ) 
                    
                    DO k = 1, num_dim
                       
                       IF ( k == i ) THEN
                          
                          this%x(k,j) = this%x(k,j) - length(k)
                          
                       ELSE
                          
                          this%x(k,j) = &
                               MODULO(this%x(k,j) - shear_length(k,2*i), &
                               length(k) )
                          
                          this%v(k,j) = this%v(k,j) + &
                               ( shear_v(k,2*i-1) - shear_v(k,2*i) )
                          
                       END IF ! k == i
                       
                    END DO ! k = 1, num_dim
                    
                    EXIT ! Found one dimension, quit.
                    
                 END SELECT ! bcdef(2*i)
                 
              END IF ! x >= max_phys
              
           END DO ! i : num_dim
           
        END DO ! j : num_colloid
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(length)) THEN
           DEALLOCATE(length) 
        END IF
        
        IF(ASSOCIATED(shear_length)) THEN
           DEALLOCATE(shear_length) 
        END IF
        
        IF(ASSOCIATED(shear_v)) THEN
           DEALLOCATE(shear_v) 
        END IF
      
        
        RETURN
        
      END SUBROUTINE colloid_adjust_colloid
      
      
