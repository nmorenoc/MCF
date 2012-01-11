      SUBROUTINE colloid_compute_accumulation_matrix(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_accumulation_matrix
        !----------------------------------------------------
        !
        ! Purpose     : Using current rotaiton matrix to compute
        !               accumulative rotation matrix for colloid.
        !
        ! Referecen   : Chen et al. 
        !               Physics of Fluids, 18, 103605, 2006.
        !
        ! Revision    : V.01  2.12.2011
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
        
        INTEGER                         :: dim,i
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        dim           = this%num_dim        
        
        IF ( this%rotate ) THEN
           
           DO i = 1, this%num_colloid
           
              this%acc_matrix(:,:,i) = &
                   MATMUL(this%rot_matrix(:,:,i), &
                   this%acc_matrix(:,:,i) )
              
              ! calculate accumulative rotation vector after acc_matrix
              
           END DO ! i = 1, num_colloid
           
        END IF ! rotate
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_compute_accumulation_matrix
      
