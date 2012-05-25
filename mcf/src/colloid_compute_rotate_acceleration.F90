      SUBROUTINE colloid_compute_rotate_acceleration(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_rotate_acceleration
        !----------------------------------------------------
        !
        ! Purpose     : Compute colloid accelerations of
        !               rotation.
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
        

        this%alpha(1:3,1:num,1:itype) = 0.0_MK
        
#else
        
        !----------------------------------------------------
        ! Save the accelerations at previous time steps.
        !----------------------------------------------------
        
        i = itype 
        
        DO WHILE ( i >= 2 )
        
           this%alpha(1:3,1:num,i) = this%alpha(1:3,1:num,i-1)
           
           i = i -1
           
        END DO
        
        !----------------------------------------------------
        ! Calculate current rotating accelerations.
        !----------------------------------------------------
    
        IF( this%rotate ) THEN
           
           IF ( dim == 2 ) THEN
              
              !----------------------------------------------
              ! For 2D rotation, no difference for different
              ! shapes, as they all move with z-xis.
              !----------------------------------------------              
              
              this%alpha(3,1:num,1) = &
                   this%torque(3,1:num) / this%mmi(3,1:num)
              
           ELSE IF ( dim == 3 )  THEN
              
              !----------------------------------------------
              ! For dD rotation, differences for different
              ! shapes have to be taken care, 
              ! as momentum of inertia is different for
              ! different shapes at different orientations.
              !----------------------------------------------              
              
              DO i = 1, num
                 
                 SELECT CASE( this%shape(i) )
                    
                 CASE ( mcf_colloid_shape_sphere )
                    
                    !----------------------------------------
                    ! For sphere, moment of inertia tensor
                    ! is diagonal with same element,
                    ! no need to transfer.
                    !----------------------------------------
                    
                    this%alpha(1:3,1:num,1) = &
                         this%torque(1:3,1:num) / this%mmi(1:3,1:num)
                    
                 CASE ( mcf_colloid_shape_ellipsoid, &
                      mcf_colloid_shape_dicolloid )
                    
                    !----------------------------------------
                    ! transfer the angular velocity and torque
                    ! into body attached frame oxyz.
                    !----------------------------------------
                    
                    omega1(1:3)=&
                         MATMUL(this%acc_matrix(1:3,1:3,i),&
                         this%omega(1:3,i,1))
                    torque1(1:3)=&
                         MATMUL(this%acc_matrix(1:3,1:3,i),&
                         this%torque(1:3,i))
                    alpha1(1) =  ( torque1(1) - &
                         (this%mmi(3,i)-this%mmi(2,i))*omega1(2)*omega1(3))/&
                         this%mmi(1,i)
                    alpha1(2) =  ( torque1(2) - &
                         (this%mmi(1,i)-this%mmi(3,i))*omega1(1)*omega1(3))/&
                         this%mmi(2,i)
                    alpha1(3) =  ( torque1(3) - &
                         (this%mmi(2,i)-this%mmi(1,i))*omega1(1)*omega1(2))/&
                         this%mmi(3,i)

                    this%alpha(1:3,i,1) = &
                         MATMUL(TRANSPOSE(this%acc_matrix(1:3,1:3,i)),&
                         alpha1(1:3))
                    
                 CASE ( mcf_colloid_shape_star )
                    
                    PRINT *, "colloid_compute_accelerate: ", &
                         "star shape is not available in 3D!"
                    stat_info = -1
                    GOTO 9999
                    
                 CASE  DEFAULT
                    
                    PRINT *, "colloid_compute_acceleration: ", &
                         "No such shape in 3D !"
                    stat_info = -1
                    GOTO 9999
                    
                 END SELECT ! shape
                 
              END DO ! i = 1, num
              
           END IF ! dim = 2
           
        ELSE 
           
           this%alpha(1:3,1:num,1) = 0.0_MK
           
        END IF ! rotate
      
    
        
        !----------------------------------------------------
        ! In context of symmetry boundaries, we have to 
        ! eliminate the force on the direction normal to 
        ! symmetry boundaries and torque.
        !----------------------------------------------------
        
        Do i = 1, dim
           
           IF ( this%bcdef(2*i-1) == ppm_param_bcdef_symmetry .AND.&
                this%bcdef(2*i) == ppm_param_bcdef_symmetry ) THEN
              
              this%alpha(1:3,1:num,1) = 0.0_MK              
              
           END IF
           
        END DO ! i = 1, dim
        
#endif        
        
9999    CONTINUE
        
        
        RETURN          
        
      END SUBROUTINE colloid_compute_rotate_acceleration
      
      
      
