      SUBROUTINE colloid_compute_rotation_matrix(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_rotation_matrix
        !----------------------------------------------------
        !
        ! Purpose     : Using current rotaiton vector to 
        !               compute current rotation matrix 
        !               for colloids at this time step.
        !
        ! Referecen   : Chen et al. 
        !               Physics of Fluids, 18, 103605, 2006.
        !
        ! Revision    : V.01  23.08.2010
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
        
        TYPE(Colloid), INTENT(OUT)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: dim,i
        REAL(MK), DIMENSION(3,3)        :: rot_matrix
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0        
        dim           = this%num_dim
        
        IF ( this%rotate ) THEN
           
           DO i = 1, this%num_colloid
              
              rot_matrix(:,:) = 0.0_MK
              
              CALL tool_rotation_matrix(this%tool,&
                   dim,this%rot_vector(1:3,i),this%rot_vector(4,i),&
                   rot_matrix(1:3,1:3),stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "colloid_compute_rotation_matrix: ", &
                      "Using tool_rotation_matrix failed! "
                 stat_info = -1
                 GOTO 9999                 
              END IF
              
              this%rot_matrix(:,:,i) = rot_matrix(:,:)
              
              !this%acc_matrix(1:dim,1:dim,i) = &
              !     MATMUL(rot_matrix(1:dim,1:dim), &
              !     this%acc_matrix(1:dim,1:dim,i) )
                   
              ! calculate accumulative rotation vector after acc_matrix
              
           END DO ! i = 1, num_colloid
           
        END IF ! rotate
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_compute_rotation_matrix
      
