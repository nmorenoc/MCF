      SUBROUTINE colloid_cartesian_ellipsoid_shortestD(&
           this,a,b,c,p,x,d,stat_info)
        !----------------------------------------------------
        ! Find the shortest distance from point 
        ! A(p(1:3)) to the surface of an ellipsoid
        !
        ! Input :
        ! center of ellipsoid is assumed to be at origin (0,0,0).
        !
        ! a   : semimajor axis along x coordinate.
        ! b   : semiminor axis along y coordinate.
        ! c   : semininor axis along z coordinate.
        !        a>=b>=c must be satisfied.
        !
        ! A is a point which can be either inside or 
        !        outside of the ellipse.
        !
        ! Output : 
        !(x(1:3):  coordinate for a point B on ellipsoid which 
        !          has shortest distance to A.
        ! d      : shortest distance from A to B.
        !
        !
        ! Remark :           
        !          1)use the symmetry, 
        !          flip A to the first Octant.
        !
        !          2) we consider for the moment only prolate
        !             and oblate.
        !----------------------------------------------------
        ! Reference :
        ! <<Distance from a Point to an Ellipse in 2D>>
        ! David Eberly,
        ! Geometric Tools, LLC, 2008.
        ! http://www.geometrictools.com/
        !----------------------------------------------------
        
        TYPE(colloid), INTENT(IN)               :: this
        REAL(MK),INTENT(IN)                     :: a
        REAL(MK),INTENT(IN)                     :: b
        REAL(MK),INTENT(IN)                     :: c
        REAL(MK),DIMENSION(:),INTENT(IN)        :: p
        REAL(MK),DIMENSION(:),INTENT(OUT)       :: x
        REAL(MK), INTENT(OUT)                   :: d
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here.
        !
        ! for prolate,
        ! point A and axis x form a plane,
        ! which cut prolate with an ellipse.
        !
        ! phi: angle between ellipse and x-y plane.
        !----------------------------------------------------
        
        INTEGER                         :: i
        LOGICAL,DIMENSION(3)            :: flip        
        REAL(MK),DIMENSION(3)           :: u
        REAL(MK)                        :: phi
        REAL(MK),DIMENSION(3)           :: ax
        REAL(MK),DIMENSION(3,3)         :: rot_m
        
        REAL(MK),DIMENSION(2)           :: pseu_p
        REAL(MK),DIMENSION(2)           :: pseu_x

        
        INTEGER                         :: stat_info_sub
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        flip(:) = .FALSE.
        
        !----------------------------------------------------
        ! Flip point to non-negative by symmetry,
        ! (first octant).
        !----------------------------------------------------
        
        u(1:3) = p(1:3)
        
        DO i = 1, 3
           
           IF ( u(i) < 0.0_MK ) THEN
              
              u(i) = -u(i)
              flip(i) = .TRUE.
              
           END IF
           
        END DO
        
        
        !----------------------------------------------------
        ! oblate
        !----------------------------------------------------
        
        IF ( ABS(a-b) < mcf_machine_zero .AND. &
             c < a ) THEN
         
           phi     = colloid_polar_angle(u(1),u(2))
           ax(1:2) = 0.0_MK
           ax(3)   = 1.0_MK
           
           CALL tool_rotation_matrix(this%tool,&
                3,ax,phi,rot_m,stat_info_sub)
           
           IF ( stat_info_sub /=  0 ) THEN
              PRINT *, __FILE__, ":", __LINE__
              stat_info = -1
              GOTO 9999
           END IF
           
           pseu_p(1) = SQRT(u(1)**2+u(2)**2)
           pseu_p(2) = u(3)
           
           CALL colloid_cartesian_ellipse_shortestD(&
                a,c,0.0_MK,pseu_p(1),pseu_p(2),&
                pseu_x(1),pseu_x(2),d,stat_info_sub)
           
           IF ( stat_info_sub /=  0 ) THEN
              PRINT *, __FILE__, ":", __LINE__
              stat_info = -1
              GOTO 9999
           END IF

           u(1) = pseu_x(1)
           u(2) = 0.0_MK
           u(3) = pseu_x(2)
           
           
           x(1:3) = MATMUL(rot_m(1:3,1:3),u(1:3))
           
           !-------------------------------------------------
           ! prolate
           !-------------------------------------------------
           
        ELSE  IF ( ABS(b-c) < mcf_machine_zero .AND. &
             b < a ) THEN
           
           phi     = colloid_polar_angle(u(2),u(3))
           ax(1)   = 1.0_MK
           ax(2:3) = 0.0_MK
           
           CALL tool_rotation_matrix(this%tool,&
                3,ax,phi,rot_m,stat_info_sub)
           
           IF ( stat_info_sub /=  0 ) THEN
              PRINT *, __FILE__, ":", __LINE__
              stat_info = -1
              GOTO 9999
           END IF
           
           pseu_p(1) = u(1)
           pseu_p(2) = SQRT(u(2)**2+u(3)**2)
           
           CALL colloid_cartesian_ellipse_shortestD(&
                a,b,0.0_MK,pseu_p(1),pseu_p(2),&
                pseu_x(1),pseu_x(2),d,stat_info_sub)
           
           IF ( stat_info_sub /=  0 ) THEN
              PRINT *, __FILE__, ":", __LINE__
              stat_info = -1
              GOTO 9999
           END IF
           
           u(1:2) = pseu_x(1:2)
           u(3)   = 0.0_MK
           
           x(1:3) = MATMUL(rot_m(1:3,1:3),u(1:3))
           
        ELSE
           
           PRINT *, __FILE__, ":", __LINE__, "no such shape!"
           stat_info = -1
           GOTO 9999
           
        END IF
        
        !----------------------------------------------------
        ! Flip back the result
        !----------------------------------------------------
        DO i = 1, 3
           
           IF ( flip(i) ) THEN
              
              x(i) = -x(i)
              
           END IF
           
        END DO
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  colloid_cartesian_ellipsoid_shortestD
      
      
    
