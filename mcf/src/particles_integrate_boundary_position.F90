      SUBROUTINE particles_integrate_boundary_position(this,&
           dt,accuracy_order,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_integrate_boundary_position
        !----------------------------------------------------
        !
        ! Purpose     : Integrate position of boundary 
        !               particles for boundary conditon
        !               with required accuracy.
        !
        ! Reference   :
        !
        ! Remark      : 
        !
        ! Revision    : V0.1  22.11.2010, original version.
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
        !  Arguments
        !
        !  this           : an object of Particles Class.
        !  num            : number of particles needed to be updated,
        !                   i.e. first num particles in this%x 
        !                   are operated.
        !  dt             : time step.
        !  accuracy_order : accuracy required.
        !  stat_info      : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        REAL(MK),INTENT(IN)                     :: dt
        INTEGER, INTENT(IN)                     :: accuracy_order
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------

        INTEGER                                 :: dim,i,ip
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        
        dim       = this%num_dim
        
#ifdef __PARTICLES_POSITION_FIXED
#else

#ifdef __WALL_FIXED
#else
        !----------------------------------------------------
        ! Select different accuracy oder
        ! 1 only velocity contribution.
        ! 2 velocity + acceleration contribution.
        !----------------------------------------------------

        
        SELECT CASE (accuracy_order)
           
        CASE (1)
           
           DO i =1, this%num_part_wall_solid_real

              ip  = this%part_wall_solid_real_list(1,i)
              
              this%x(1:dim,ip) = &
                   this%x(1:dim,ip) + &
                   this%v(1:dim,ip) * dt
                 
           END DO
           
        CASE (2)
           
           DO i = 1, this%num_part_wall_solid_real
              
              this%x(1:dim,ip) = &
                   this%x(1:dim,ip) + &
                   this%v(1:dim,ip) * dt
              ! + &
              !     0.5_MK * this%f(1:dim,ip) * dt**2
                 
           END DO
           
           !-------------------------------------------------
           ! Other integration not available.
           !-------------------------------------------------
           
        CASE DEFAULT
           PRINT *, "particles_integrate_boundary_position : ", &
                "Order of accuracy not available !"
           stat_info = -1
           GOTO 9999
           
        END SELECT ! accuracy_order
#endif

#endif
        
9999    CONTINUE      
        
        RETURN
        
      END SUBROUTINE particles_integrate_boundary_position
      
