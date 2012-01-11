      SUBROUTINE particles_compute_pressure(this,num,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_compute_pressure
        !----------------------------------------------------
        !
        ! Purpose     : Computing pressure of particles.
        !                 
        !      
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.2 08.07 2009,
        !               check again the work flow is correct
        !               and supply with more comments.
        !
        !               V0.1 01.04 2009, original version.
        !-----------------------------------------------------
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
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
      	! Local variables starts here :
      	!----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: rhs_density_type
        REAL(MK), DIMENSION(:), POINTER         :: t_p
        
        
        !----------------------------------------------------
      	! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        NULLIFY(t_p)
        
        
        IF( num > this%num_part_all ) THEN
           PRINT *, "particles_compute_pressure : ", &
                "num > num_part_all, wrong !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        rhs_density_type = &
             control_get_rhs_density_type(this%ctrl,stat_info_sub)
        
             
        SELECT  CASE( rhs_density_type )
           
        CASE (1)
           
           CALL stateEquation_compute_pressure(this%stateEquation,&
                this%rho(1:num),t_p,num,stat_info_sub)
           
        CASE (2)
           
           CALL stateEquation_compute_pressure(this%stateEquation,&
                this%rho(1:num)*this%m(1:num),&
                t_p,num,stat_info_sub)
           
        END SELECT
        
        IF( stat_info_sub /=0 ) THEN
           PRINT *, "particles_compute_pressure : ", &
                "stateEquation_compute_pressure has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        IF(ASSOCIATED(this%p)) THEN
           DEALLOCATE(this%p)
        END IF
        
        ALLOCATE(this%p(1:num))
        
        this%p(1:num) = t_p(1:num)
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(t_p)) THEN
           DEALLOCATE(t_p)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_compute_pressure
      
