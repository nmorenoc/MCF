      SUBROUTINE colloid_cartesian_ellipse_shortestD(&
           a,b,phi,p,q,x,y,d,stat_info)
        !----------------------------------------------------
        ! Find the shortest distance from point A(p,q) to
        ! the surface of an ellipse.
        !
        ! Input :
        ! center of ellipse is assumed to be at origin (0,0).
        !
        ! a   : semimajor axis along x coordinate.
        ! b   : semiminor axis along y coordinate.
        !       a>=b must be satisfied.
        ! phi : ellipse rotate phi radian to x direction.
        !       counter-clockwise is positive.
        ! A(p,q) is a point which can be either inside or 
        !        outside of the ellipse.
        !
        ! Output : 
        !(x,y): coordinate for a point B on ellipse which 
        !       has shortest distance to A(u,v)
        ! d   : shortest distance from A to B.
        !
        !
        ! Remark : 1)rotate point A(p,q) phi radian clockwise,
        !          
        !          2)use the symmetry of ellipse, assume
        !          flip A(x,v) to the first quadrant.
        !
        !          3)find shortest distance from A(u,v) to
        !          the ellipse x**2/a**2 + y**2/b**2 = 1 at 
        !          point B(x',v') on the surface.
        !
        !          4)flip B(x',v') back using symmetry.
        !
        !          5)rotate point B(x',v') counter-clockwise
        !          phi radian to get B(x,v).
        !
        !----------------------------------------------------
        ! Reference :
        ! <<Distance from a Point to an Ellipse in 2D>>
        ! David Eberly,
        ! Geometric Tools, LLC, 2008.
        ! http://www.geometrictools.com/
        !----------------------------------------------------
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: phi        
        REAL(MK), INTENT(IN)            :: p
        REAL(MK), INTENT(IN)            :: q
        REAL(MK), INTENT(OUT)           :: x
        REAL(MK), INTENT(OUT)           :: y
        REAL(MK), INTENT(OUT)           :: d
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here.
        !----------------------------------------------------
        
        LOGICAL                         :: flip_x
        LOGICAL                         :: flip_y
        
        REAL(MK)                        :: theta
        REAL(MK)                        :: u
        REAL(MK)                        :: v        
        REAL(MK)                        :: len

        REAL(MK)                        :: t0
        REAL(MK)                        :: t1
        REAL(MK)                        :: f, df
        INTEGER                         :: iter
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        flip_x = .FALSE.
        flip_y = .FALSE.
        
        !----------------------------------------------------
        ! Rotate A(p,q) to A(u,v) by phi radian clockwise
        !----------------------------------------------------
        
        theta = colloid_polar_angle(p,q)
        len   = SQRT(p**2+q**2)
        
        u = len*COS(theta-phi)
        v = len*SIN(theta-phi)
        
        !----------------------------------------------------
        ! Flip point A(u,v) to non-negative u,v by symmetry,
        ! i.e.,(0,0),(u,0),(0,v) or (u,v)(first quadrant).
        !----------------------------------------------------
        
        IF ( u < 0.0_MK ) THEN
           
           u      = -u
           flip_x = .TRUE.
           
        END IF
        
        IF ( v < 0.0_MK ) THEN
           
           v      = -v
           flip_y = .TRUE.
           
        END IF
        
        !----------------------------------------------------
        ! Check where is point A(u,v).
        !----------------------------------------------------
        
        IF ( ABS(u) < mcf_machine_zero) THEN
           
           !-------------------------------------------------
           ! point A(u,v) is at origin (0,0) or
           ! point A(u,v) is along Y axis.           
           !-------------------------------------------------
           
           x = 0.0_MK
           y = b
           d = ABS(b-v)
           
        ELSE IF ( ABS(v) < mcf_machine_zero ) THEN
        
           !-------------------------------------------------
           ! point A(u,v) is along X axis.
           !-------------------------------------------------
           
           IF ( u < a-b**2/a ) THEN
              
              x = a**2*u/(a**2-b**2)
              y = b*SQRT(1-(x/a)**2)
              d = b*SQRT(1-u**2/(a**2-b**2))
              
           ELSE
              
              x = a
              y = 0.0_MK
              d = ABS(u-a)
              
           END IF
           
        ELSE
           
           !-------------------------------------------------
           ! point A(u,v) is in first quadrant
           ! (non-zero cooridnate).
           !-------------------------------------------------
           
           !-------------------------------------------------
           ! Since the value of f(t) is mono-increasing,
           ! therefore, we can use  Newton's Method for
           ! finding the root.
           ! t0 is the first guess.
           !-------------------------------------------------
           
           t0 = b*v-b**2
           f  = colloid_cartesian_ellipse_F(a,b,u,v,t0)
           
           iter = 0
           
           DO WHILE( ABS(f) > mcf_machine_zero .AND. &
                iter < 20 )
              
              df = colloid_cartesian_ellipse_dF(a,b,u,v,t0)
              
              t1 = t0 - f/df
              
              f = colloid_cartesian_ellipse_F(a,b,u,v,t1)
              
              t0 = t1
              iter = iter + 1
              
           END DO
           
           x = a**2*u/(t0+a**2)
           y = b**2*v/(t0+b**2)
           d = SQRT((x-u)**2+(v-y)**2)
           
        END IF
        
        !----------------------------------------------------
        ! Flip the point on ellipse which has shortest
        ! distance from A(u,v) using symmetry.
        !----------------------------------------------------
        
        IF ( flip_x ) THEN
           x = -x
        END IF
        
        IF ( flip_y ) THEN
           y = -y
        END IF
        
        !----------------------------------------------------
        ! Rotate counter-closewise (x,y) by phi radian.
        !----------------------------------------------------
        
        theta = colloid_polar_angle(x,y)
        len   = SQRT(x**2+y**2)
        x     = len*COS(theta+phi)
        y     = len*SIN(theta+phi)
        
        RETURN
        
      END SUBROUTINE  colloid_cartesian_ellipse_shortestD


      REAL(MK) FUNCTION colloid_cartesian_ellipse_F(a,b,u,v,t)
        !----------------------------------------------------
        ! F(t) = (a*u/(t+a**2))**2 + (b*v/(t+b**2))**2 - 1
        !----------------------------------------------------

        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        REAL(MK), INTENT(IN)            :: t
        
        colloid_cartesian_ellipse_F = &
             (a*u/(t+a**2))**2 + &
             (b*v/(t+b**2))**2 - 1
        
      END FUNCTION colloid_cartesian_ellipse_F
      
      
      REAL(MK) FUNCTION colloid_cartesian_ellipse_dF(a,b,u,v,t)
        !----------------------------------------------------
        ! First derivative of F(t)
        ! F1 = -2*a**2*u**2/(t+a**2)**3-2*b**2*v**2/(t+b**2)**3
        !----------------------------------------------------

        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        REAL(MK), INTENT(IN)            :: t
        
        
        colloid_cartesian_ellipse_dF = &
             -2.0_MK*a**2*u**2 / (t+a**2)**3 - &
             2.0_MK*b**2*v**2 / (t+b**2)**3
        
        
      END FUNCTION colloid_cartesian_ellipse_dF
      
      
      REAL(MK) FUNCTION colloid_cartesian_ellipse_ddF(a,b,u,v,t)
        !----------------------------------------------------
        ! Second derivative of F(t)
        ! F2 = 6*a**2*u**2/(t+a**2)**4+6*b**2*v**2/(t+b**2)**4
        !----------------------------------------------------

        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        REAL(MK), INTENT(IN)            :: t
        
        
        colloid_cartesian_ellipse_ddF = &
             6.0_MK*a**2*u**2 / (t+a**2)**4 + &
             6.0_MK*b**2*v**2 / (t+b**2)**4 
        
        
      END FUNCTION colloid_cartesian_ellipse_ddF
      
