      SUBROUTINE colloid_finalize(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_finalize
        !----------------------------------------------------
        !
        ! Purpose     : Destrutor of Class Colloid.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 01.03.2009, original version.
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
        
        TYPE(Colloid),POINTER                 :: this
        INTEGER,INTENT(OUT)                   :: stat_info
        
        
        stat_info = 0
        
        IF( ASSOCIATED(this%body_force) ) THEN
           DEALLOCATE(this%body_force)             
        END IF
        
        IF( ASSOCIATED(this%shape) ) THEN
           DEALLOCATE(this%shape)
        END IF
        
        IF( ASSOCIATED(this%radius) ) THEN
           DEALLOCATE(this%radius)
        END IF

        IF( ASSOCIATED(this%freq) ) THEN
           DEALLOCATE(this%freq)
        END IF
        
        IF( ASSOCIATED(this%m) ) THEN
           DEALLOCATE(this%m)             
        END IF
        
        IF( ASSOCIATED(this%mmi) ) THEN
           DEALLOCATE(this%mmi)
        END IF
        
        IF( ASSOCIATED(this%x) ) THEN
           DEALLOCATE(this%x)
        END IF

        IF( ASSOCIATED(this%v) ) THEN
           DEALLOCATE(this%v)
        END IF
        
#if __DRAG_PART
        IF(ASSOCIATED(this%drag_lub) ) THEN
           DEALLOCATE(this%drag_lub)
        END IF
        
        IF(ASSOCIATED(this%drag_repul) ) THEN
           DEALLOCATE(this%drag_repul)
        END IF
#endif

        IF(ASSOCIATED(this%drag) ) THEN
           DEALLOCATE(this%drag)
        END IF
        
        IF( ASSOCIATED(this%rot_vector) ) THEN
           DEALLOCATE(this%rot_vector)
        END IF

        IF( ASSOCIATED(this%acc_vector) ) THEN
           DEALLOCATE(this%acc_vector)
        END IF

        IF( ASSOCIATED(this%rot_matrix) ) THEN
           DEALLOCATE(this%rot_matrix)
        END IF

        IF( ASSOCIATED(this%acc_matrix) ) THEN
           DEALLOCATE(this%acc_matrix)
        END IF
        
        IF( ASSOCIATED(this%theta) ) THEN
           DEALLOCATE(this%theta)
        END IF
        
        IF( ASSOCIATED(this%omega) ) THEN
           DEALLOCATE(this%omega)
        END IF
        
        IF(ASSOCIATED(this%torque) ) THEN
           DEALLOCATE(this%torque)
        END IF
        
        IF( ASSOCIATED(this%num_physical_part) ) THEN
           DEALLOCATE(this%num_physical_part)
        END IF
        
        IF( ASSOCIATED(this%num_numerical_part) ) THEN
           DEALLOCATE(this%num_numerical_part)
        END IF
        
        IF( ASSOCIATED(this%f) ) THEN
           DEALLOCATE(this%f)             
        END IF

        IF( ASSOCIATED(this%alpha) ) THEN
           DEALLOCATE(this%alpha)
        END IF
        
        IF( ASSOCIATED(this%k_energy) ) THEN
           DEALLOCATE(this%k_energy)             
        END IF
        
        IF( ASSOCIATED(this%mom) ) THEN
           DEALLOCATE(this%mom)             
        END IF
        
        IF( ASSOCIATED(this%mom_tot) ) THEN
           DEALLOCATE(this%mom_tot)             
        END IF

        IF( ASSOCIATED(this%min_phys) ) THEN
           DEALLOCATE(this%min_phys)
        END IF
        
        IF( ASSOCIATED(this%max_phys) ) THEN
           DEALLOCATE(this%max_phys)
        END IF
     
        IF( ASSOCIATED(this%min_phys_t) ) THEN
           DEALLOCATE(this%min_phys_t)
        END IF
        
        IF( ASSOCIATED(this%max_phys_t) ) THEN
           DEALLOCATE(this%max_phys_t)
        END IF

        IF( ASSOCIATED(this%bcdef) ) THEN
           DEALLOCATE(this%bcdef)
        END IF

        IF( ASSOCIATED(this%x_image) ) THEN
           DEALLOCATE(this%x_image)
        END IF
        
        IF( ASSOCIATED(this%v_image) ) THEN
           DEALLOCATE(this%v_image)
        END IF

        RETURN          
        
      END SUBROUTINE colloid_finalize
      
