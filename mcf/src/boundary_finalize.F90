!-------------------------------------------------
! Destructors of Class Boundary
!-------------------------------------------------

      SUBROUTINE boundary_finalize(this,stat_info)

        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0
        
        IF (ASSOCIATED(this%bcdef)) THEN
           DEALLOCATE(this%bcdef)
        END IF
        
        IF (ASSOCIATED(this%shear_rate)) THEN
           DEALLOCATE(this%shear_rate)
        END IF

        IF (ASSOCIATED(this%shear_length)) THEN
           DEALLOCATE(this%shear_length)
        END IF
        
        IF (ASSOCIATED(this%shear_type)) THEN
           DEALLOCATE(this%shear_type)
        END IF
        
        IF (ASSOCIATED(this%shear_v0)) THEN
           DEALLOCATE(this%shear_v0)
        END IF
        
        IF (ASSOCIATED(this%shear_v)) THEN
           DEALLOCATE(this%shear_v)
        END IF

        IF (ASSOCIATED(this%shear_freq)) THEN
           DEALLOCATE(this%shear_freq)
        END IF

        IF (ASSOCIATED(this%drag)) THEN
           DEALLOCATE(this%drag)
        END IF
        
        IF (ASSOCIATED(this%drag_p)) THEN
           DEALLOCATE(this%drag_p)
        END IF
        IF (ASSOCIATED(this%drag_v)) THEN
           DEALLOCATE(this%drag_v)
        END IF
        IF (ASSOCIATED(this%drag_r)) THEN
           DEALLOCATE(this%drag_r)
        END IF

        
        RETURN
        
      END SUBROUTINE boundary_finalize
