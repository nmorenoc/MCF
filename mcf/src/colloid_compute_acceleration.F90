      SUBROUTINE colloid_compute_acceleration(this,time,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_acceleration
        !----------------------------------------------------
        !
        ! Purpose     : Compute colloid accelerations,
        !               i.e., translation and rotation.
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
        !               torque needs to be decomposed.
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
        REAL(MK), INTENT(IN)                    :: time
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                                 :: dim,num
        INTEGER                                 :: i
        REAL(MK), DIMENSION(3)                  :: omega1
        REAL(MK), DIMENSION(3)                  :: torque1
        REAL(MK), DIMENSION(3)                  :: alpha1

        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        
        dim   = this%num_dim
        num   = this%num_colloid
        
        
#ifdef __COLLOID_NOACCE
        
        !----------------------------------------------------
        ! No acceleration
        !        
        ! Solution to keep colloid moving/rotating
        ! in a constant velocity, i.e., no translational or
        ! rotational acceleration.
        !----------------------------------------------------
        
        this%f(1:dim,1:num)   = 0.0_MK
        this%alpha(1:3,1:num) = 0.0_MK
        
#else
        
        !----------------------------------------------------
        ! Calculate translating accelerations.
        !----------------------------------------------------
        
        IF( this%translate ) THEN
           
           DO i = 1, dim
              
              this%f(i,1:num) = &
                   this%drag(i,1:num) / this%m(1:num)
              
           END DO
           
        ELSE
           
           this%f(1:dim,1:num) = 0.0_MK
           
        END IF
        
        !----------------------------------------------------
        ! No movement in y-direction.
        !----------------------------------------------------
        
        !IF ( time < 30 ) THEN
        !   this%f(2,1:num) = 0.0_MK
        !END IF
     
        !----------------------------------------------------
        ! Calculate rotating accelerations.
        ! Note that torque and its accleration are 3D.
        !----------------------------------------------------
        
        IF( this%rotate ) THEN
           
           IF ( dim == 2 ) THEN
              
              this%alpha(3,1:num) = &
                   this%torque(3,1:num) / this%mmi(3,1:num)
              
              
           ELSE IF ( dim == 3 )  THEN
              
              DO i = 1, num
                 
                 SELECT CASE( this%shape(i) )
                    
                 CASE ( mcf_colloid_shape_sphere )
                    !----------------------------------------
                    ! For sphere, moment of inertia tensor
                    ! is diagonal with same element,
                    ! no need to transfer.
                    !----------------------------------------

                    this%alpha(1:3,1:num) = &
                      this%torque(1:3,1:num) / this%mmi(1:3,1:num)
              
                 CASE ( mcf_colloid_shape_ellipsoid, &
                      mcf_colloid_shape_dicolloid )
                    !----------------------------------------
                    ! transfer the angular velocity and torque
                    ! into body attached frame oxyz.
                    !----------------------------------------
                    
                    omega1(1:3)=&
                         MATMUL(this%acc_matrix(1:3,1:3,i),&
                         this%omega(1:3,i))
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

                    this%alpha(1:3,i) = &
                         MATMUL(TRANSPOSE(this%acc_matrix(1:3,1:3,i)),&
                         alpha1(1:3))
                    
                 CASE ( mcf_colloid_shape_star )
                    PRINT *, "colloid_compute_accelerate: ", &
                         "star shape is not available in 3D !"
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
           
           this%alpha(1:3,1:num) = 0.0_MK
           
        END IF ! rotate
        
        
        !----------------------------------------------------
        ! In context of symmetry boundaries, we have to 
        ! eliminate the force on the direction normal to 
        ! symmetry boundaries and torque.
        !----------------------------------------------------
        
        Do i = 1, dim
           
           IF ( this%bcdef(2*i-1) == ppm_param_bcdef_symmetry .AND.&
                this%bcdef(2*i) == ppm_param_bcdef_symmetry ) THEN
              
              this%f(i,1:num) = 0.0_MK
              
              this%alpha(1:3,1:num) = 0.0_MK              
              
           END IF
           
        END DO ! i = 1, dim
        
#endif        
        
        
9999    CONTINUE
        
        
        RETURN          
        
      END SUBROUTINE colloid_compute_acceleration
      
      
      
