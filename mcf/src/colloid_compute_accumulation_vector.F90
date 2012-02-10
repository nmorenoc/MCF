      SUBROUTINE colloid_compute_accumulation_vector(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_accumulation_vector
        !----------------------------------------------------
        !
        ! Purpose     : Compute the accumulative rotation vector
        !               using accumulative rotation matrix.
        !
        ! Remark      : Colloid are modeled as rigid body.
        !
        ! Reference  : Chen et. al. 2006, physics of fluids
        !              wikipedia
        !
        ! Revision   : V0.1  2.12.2011, original.
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
        REAL(MK),DIMENSION(3)           :: axis
        REAL(MK)                        :: phi
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        dim           = this%num_dim

#if __PARTICLES_POSITION_FIXED
#else        
        IF ( this%rotate ) THEN

           DO i = 1, this%num_colloid
                 
              CALL tool_rotation_vector(this%tool, &
                   dim,this%acc_matrix(1:3,1:3,i),&
                   axis(1:3),phi,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "colloid_compute_accumulation_vector ; ", &
                      "Using tool_rotation_vector failed ! "
                 stat_info = -1
                 GOTO 9999
              END IF

              this%acc_vector(1:3,i) = axis(1:3)
              this%acc_vector(4,i)   = phi
              
           END DO ! i = 1, num_colloid
           
        END IF ! rotate
        
#endif
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_compute_accumulation_vector
      
