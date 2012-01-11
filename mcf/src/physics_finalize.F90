      SUBROUTINE physics_finalize(this,stat_info)
        !----------------------------------------------------
        !  Subroutine   :  physics_finalize
        !----------------------------------------------------
        !
        !  Purpose      : Destrutor of Class Physics.
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : V0.1 01.03.2009, original version.
        !
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        TYPE(Physics),INTENT(INOUT)     :: this
        INTEGER,INTENT(OUT)             :: stat_info
        
        !------------------------
        ! Finalize the physics object
        !-----------------------
        INTEGER                         :: stat_info_sub   
        
        stat_info     = 0        
        stat_info_sub = 0        

        IF(ASSOCIATED(this%min_phys)) THEN
           DEALLOCATE(this%min_phys)
        END IF

        IF(ASSOCIATED(this%max_phys)) THEN
           DEALLOCATE(this%max_phys)
        END IF
        
        IF(ASSOCIATED(this%min_phys_t)) THEN
           DEALLOCATE(this%min_phys_t)
        END IF

        IF(ASSOCIATED(this%max_phys_t)) THEN
           DEALLOCATE(this%max_phys_t)
        END IF
      
        IF(ASSOCIATED(this%num_part_dim)) THEN
           DEALLOCATE(this%num_part_dim)
        END IF

        IF(ASSOCIATED(this%num_part_dim_t)) THEN
           DEALLOCATE(this%num_part_dim_t)
        END IF
        
        IF(ASSOCIATED(this%dx)) THEN
           DEALLOCATE(this%dx)
        END IF

        IF(ASSOCIATED(this%body_force)) THEN
           DEALLOCATE(this%body_force)
        END IF

        IF(ASSOCIATED(this%body_force_d)) THEN
           DEALLOCATE(this%body_force_d)
        END IF

        IF(ASSOCIATED(this%eval)) THEN
           DEALLOCATE(this%eval)
        END IF

        IF(ASSOCIATED(this%evec)) THEN
           DEALLOCATE(this%evec)
        END IF

        IF(ASSOCIATED(this%bcdef)) THEN
           DEALLOCATE(this%bcdef)
        END IF
        
        IF(ASSOCIATED(this%boundary)) THEN
           CALL boundary_finalize(this%boundary,stat_info_sub)
           DEALLOCATE(this%boundary)
        END IF
        
        IF(ASSOCIATED(this%colloids)) THEN
           CALL colloid_finalize(this%colloids,stat_info_sub)
           DEALLOCATE(this%colloids)
        END IF
        
        PRINT *, "physics_finalize : ", "Finished!"
        
        RETURN
        
      END SUBROUTINE physics_finalize
     
