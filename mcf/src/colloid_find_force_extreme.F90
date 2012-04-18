      SUBROUTINE colloid_find_force_extreme(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_find_force_extreme
        !----------------------------------------------------
        !
        ! Purpose     : Find the minimal and maximal
        !               force per unit mass.
        !               
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     :
        !
        ! Revisions   : V0.1 13.10.2010, original version.
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
        ! Arguments
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(INOUT)            :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim, i, num
        REAL(MK), DIMENSION(:), POINTER         :: fa
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        

        dim = this%num_dim
        num = this%num_colloid
        
        NULLIFY(fa)
        
        ALLOCATE(fa(num))
        
        DO i = 1, num
           
           fa(i) = &
                SQRT(DOT_PRODUCT(this%f(1:dim,i,1), this%f(1:dim,i,1)))
           
        END DO
        
        this%fa_min = &
             MINVAL(fa(1:num))
        this%fa_max = &
             MAXVAL(fa(1:num))
        
        
        IF ( ASSOCIATED(fa)) THEN
           DEALLOCATE(fa)
        END IF

        RETURN
        
      END SUBROUTINE colloid_find_force_extreme
      
      
      
