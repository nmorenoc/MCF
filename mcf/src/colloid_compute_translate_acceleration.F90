      SUBROUTINE colloid_compute_translate_acceleration(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_translate_acceleration
        !----------------------------------------------------
        !
        ! Purpose     : Compute colloid accelerations of
        !               translation.
        !               
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     :
        !
        ! Revisions   : V0.2 21.11.2011, include 3D
        !               ellipsoid and dicolloid, where
        !               torque needs to be decomposed,
        !               as momentum of inertia is different
        !               at different rotating orientation.
        !
        !               V0.1 19.11.2010, original version.
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(OUT)              :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                                 :: dim,num
        INTEGER                                 :: itype,i
        REAL(MK), DIMENSION(3)                  :: omega1
        REAL(MK), DIMENSION(3)                  :: torque1
        REAL(MK), DIMENSION(3)                  :: alpha1

        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        
        dim   = this%num_dim
        num   = this%num_colloid
        itype = this%integrate_AB
        
#ifdef __COLLOID_NOACCE
        
        !----------------------------------------------------
        ! No acceleration
        !        
        ! Solution to keep colloid moving/rotating
        ! in a constant velocity, i.e., no translational or
        ! rotational acceleration.
        !----------------------------------------------------
        
        this%f(1:dim,1:num,1:itype)   = 0.0_MK
        
#else
        
        !----------------------------------------------------
        ! Save the accelerations at previous time steps.
        !----------------------------------------------------
        
        i = itype 
        
        DO WHILE ( i >= 2 )
        
           this%f(1:dim,1:num,i)  =  this%f(1:dim,1:num,i-1)                
           
           i = i -1
           
        END DO
        
        !----------------------------------------------------
        ! Calculate current translating accelerations,
        ! set zero if no translation.
        !----------------------------------------------------
        
        IF( this%translate ) THEN
           
           DO i = 1, dim
              
              this%f(i,1:num,1) = &
                   this%drag(i,1:num) / this%m(1:num)
              
           END DO
           
        ELSE
           
           this%f(1:dim,1:num,1) = 0.0_MK
           
        END IF
        
        !----------------------------------------------------
        ! In context of symmetry boundaries, we have to 
        ! eliminate the force on the direction normal to 
        ! symmetry boundaries and torque.
        !----------------------------------------------------
        
        Do i = 1, dim
           
           IF ( this%bcdef(2*i-1) == ppm_param_bcdef_symmetry .AND.&
                this%bcdef(2*i) == ppm_param_bcdef_symmetry ) THEN
              
              this%f(i,1:num,1) = 0.0_MK
              
           END IF
           
        END DO ! i = 1, dim
        
#endif        
        
9999    CONTINUE
        
        
        RETURN          
        
      END SUBROUTINE colloid_compute_translate_acceleration
      
      
      
