      MODULE COLLOID

        IMPLICIT NONE
        
        INTEGER, PARAMETER              :: MK=KIND(1.0D0)
        REAL(MK),PARAMETER              :: mcf_pi = 3.141592653589793238_MK
        REAL(MK),PARAMETER     :: mcf_machine_zero = EPSILON(mcf_pi) * 2.0_MK

      CONTAINS
        
        REAL(MK) FUNCTION angle(x,y)
          
          !--------------------------------------------------
          ! Get angle between point (x,y) and x+ direction.
          ! ACOS returns a value in [0:pi], therefore,
          ! if y coordinate is negative,
          ! we do 2*pi-angle.
          !--------------------------------------------------
          
          REAL(MK),INTENT(IN)           :: x
          REAL(MK),INTENT(IN)           :: y
          
          
          angle = ACOS(x/SQRT(x**2+y**2))
          
          IF (y < 0.0_MK) THEN
             
             angle = 2.0_MK*mcf_pi - angle
             
          END IF
          
        END FUNCTION angle
        
        
        REAL(MK) FUNCTION polar_ellipse_r(a,b,theta,phi)
          
          REAL(MK), INTENT(IN)            :: a
          REAL(MK), INTENT(IN)            :: b
          REAL(MK), INTENT(IN)            :: theta
          REAL(MK), INTENT(IN)            :: phi
          
          polar_ellipse_r = SQRT(2.0_MK) * a * b / &
               SQRT( ( b**2-a**2 )* &
               COS(2.0_MK*(theta-phi)) + &
               a**2 + b**2 )
          
        END FUNCTION polar_ellipse_r
      
      
      REAL(MK) FUNCTION polar_ellipse_r1(a,b,theta,phi)
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        polar_ellipse_r1 = &
             ( polar_ellipse_r(a,b,theta,phi) ) ** 3.0 * &
             ( b**2-a**2 ) * SIN(2.0_MK*(theta-phi)) / &
             (2.0_MK*a**2*b**2)
        
        
      END FUNCTION polar_ellipse_r1
      
      
      REAL(MK) FUNCTION polar_ellipse_r2(a,b,theta,phi)
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        REAL(MK)                        :: r
        
        r =  polar_ellipse_r(a,b,theta,phi) 
        
        polar_ellipse_r2 = 3.0_MK * &
             ( polar_ellipse_r1(a,b,theta,phi) ) ** 2 / &
             r + r**3 * ( b**2-a**2 ) * COS(2.0_MK*(theta-phi)) / &
             a**2/b**2
        
        
      END FUNCTION polar_ellipse_r2
      
      
      REAL(MK) FUNCTION polar_star_r(a,b,c,theta,phi)
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        polar_star_r = a + b * COS(c*(theta-phi))               
        
      END FUNCTION polar_star_r
      
      
      REAL(MK) FUNCTION polar_star_r1(b,c,theta,phi)
        
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        polar_star_r1 =  -b * c * SIN(c*(theta-phi))
        
      END FUNCTION polar_star_r1
      
      
      REAL(MK) FUNCTION polar_star_r2(b,c,theta,phi)
        
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        polar_star_r2 =  -b * c * c * COS(c*(theta-phi))
        
      END FUNCTION polar_star_r2
      
      
      REAL(MK) FUNCTION polar_star_F(a,b,c,theta,u,v)
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        
        REAL(MK)                        :: r
        REAL(MK)                        :: r1
        
        r  = polar_star_r(a,b,c,theta,0.0_MK)
        r1 = polar_star_r1(b,c,theta,0.0_MK)
        
        polar_star_F = r*r1 + &
             r  * ( u*SIN(theta) - v*COS(theta)) - &
             r1 * ( v*SIN(theta) + u*COS(theta)) 
        
      END FUNCTION polar_star_F
      
      
      REAL(MK) FUNCTION polar_star_F1(a,b,c,theta,u,v)
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        
        REAL(MK)                        :: r
        REAL(MK)                        :: r1
        REAL(MK)                        :: r2
        
        
        r  = polar_star_r(a,b,c,theta,0.0_MK)
        r1 = polar_star_r1(b,c,theta,0.0_MK)
        r2 = polar_star_r2(b,c,theta,0.0_MK)
        
        polar_star_F1 = r*r2 + r1**2 + &
             2.0_MK*r1*( u*SIN(theta) - v*COS(theta)) + &
             ( r-r2) *( v*SIN(theta) + u*COS(theta))
        
      END FUNCTION polar_star_F1
      
      
      REAL(MK) FUNCTION polar_star_F_root(a,b,c,u,v,&
           theta_min,theta_max,stat_info)
        
        !----------------------------------------------------
        ! Find the root of function F(theta) = 0.0
        ! in (theta_min, theta_max).
        !
        ! Procedure :
        ! We use bisection.
        !----------------------------------------------------
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        REAL(MK), INTENT(in)            :: theta_min
        REAL(MK), INTENT(in)            :: theta_max
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here.
        !
        !----------------------------------------------------
   
        REAL(MK)                        :: f1
        REAL(MK)                        :: f2
        REAL(MK)                        :: fh

        REAL(MK)                        :: theta1
        REAL(MK)                        :: theta2
        REAL(MK)                        :: thetah
        
        INTEGER                         :: iter_max
        INTEGER                         :: iter
        
        stat_info = 0
        iter_max = 50
        iter     = 0

        
        theta1 = theta_min
        theta2 = theta_max
        
        f1 = polar_star_F(a,b,c,theta1,u,v)
        f2 = polar_star_F(a,b,c,theta2,u,v)
        
        IF ( f1*f2 >= 0.0_MK ) THEN
           !PRINT *, "colloid_find_F_root : ", &
           !     "Root must be bracketed"
           !PRINT *, "u, v : ", u,v
           stat_info = -1
           GOTO 9999
        END IF
        
        thetah = (theta1 + theta2) / 2.0_MK
        fh     = polar_star_F(a,b,c,thetah,u,v)
        
        DO WHILE( ABS(fh) > mcf_machine_zero .AND. &
             iter < iter_max )
           
           iter = iter + 1
           
           IF ( fh * f1 > 0.0_MK ) THEN
              
              theta1 = thetah
              f1     = polar_star_F(a,b,c,theta1,u,v)
              
           ELSE

              theta2 = thetah
              f2     = polar_star_F(a,b,c,theta1,u,v)

           END IF
           
           thetah = (theta1 + theta2) / 2.0_MK
           fh = polar_star_F(a,b,c,thetah,u,v)
           
        END DO
        
        IF ( iter >= iter_max ) THEN
           
           polar_star_F_root = -1.0_MK
           
        ELSE
           
           polar_star_F_root = thetah
           
        END IF
        
        !PRINT *, iter
        
9999    CONTINUE
        
        RETURN
        
      END FUNCTION polar_star_F_root
      
      
      SUBROUTINE polar_star_shortestD(a,b,c,phi,p,q,x,y,d)
        
        !----------------------------------------------------
        ! Find the shortest distance between a point (p,q)
        ! to the curve given in polar coordinate 
        ! r(t) = a+b*cos(c*(t-phi));
        !
        ! Procedure :
        ! We first rotate axis by phi radian counter-clockwise,
        ! then always map point(p,q) by symmetry axis to
        ! the t=[0:pi/c] and do calculation there
        ! to find x,y,d by using polar_star_F_root
        ! (bisection).
        !
        ! Afterwards and map (x,y) back and 
        ! rotate back by phi radian.
        !----------------------------------------------------
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: phi
        REAL(MK), INTENT(IN)            :: p
        REAL(MK), INTENT(IN)            :: q
        REAL(MK), INTENT(OUT)           :: x
        REAL(MK), INTENT(OUT)           :: y
        REAL(MK), INTENT(OUT)           :: d
        
        !----------------------------------------------------
        ! Local variables start here.
        !
        ! angle_sec : angle of each section of symmetry.
        !
        !----------------------------------------------------
        
        REAL(MK)                        :: theta
        REAL(MK)                        :: u
        REAL(MK)                        :: v        
        REAL(MK)                        :: r
        REAL(MK)                        :: angle_sec
        REAL(MK)                        :: angle_min
        REAL(MK)                        :: angle_max
        INTEGER                         :: sec
        LOGICAL                         :: flip
        REAL(MK)                        :: f1,f2
        REAL(MK)                        :: r1,r2
        REAL(MK)                        :: x1,y1
        REAL(MK)                        :: x2,y2
        REAL(MK)                        :: d1,d2
        REAL(MK)                        :: theta1,theta2

        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        flip = .FALSE.
        stat_info_sub = 0
        
        !----------------------------------------------------
        ! Rotate axis by phi radian couter-clockwise.
        !----------------------------------------------------
        
        theta = angle(p,q)
        r = SQRT(p**2+q**2)
        
        theta = theta - phi
        
        !----------------------------------------------------
        ! Rotate theta to first section of symmetry.
        !----------------------------------------------------
        
        angle_sec = 2.0_MK * mcf_pi / c
        sec       = INT(theta/angle_sec)
        theta     = theta - angle_sec * sec
        
        !----------------------------------------------------
        ! Flip angle in first section due to symmetry.
        !----------------------------------------------------
        
        IF ( theta > angle_sec / 2.0_MK ) THEN
           
           flip = .TRUE.
           theta = angle_sec - theta
           
        END IF
        
        angle_min = 0.0_MK
        angle_max = mcf_pi / c
        
        !----------------------------------------------------
        ! find point(u,v) in new coordinate for (p,q).
        !----------------------------------------------------
        
        u = r*COS(theta)
        v = r*SIN(theta)
        
        !PRINT *, u,v
        
        !----------------------------------------------------
        ! check if angle_min or angle_max is the root of
        ! function F.
        !----------------------------------------------------
        
        f1 = polar_star_F(a,b,c,angle_min,u,v)
        f2 = polar_star_F(a,b,c,angle_max,u,v)
        
        IF ( ABS(f1) < mcf_machine_zero .AND. &
             ABS(f2) < mcf_machine_zero ) THEN
           
           r1 = polar_star_r(a,b,c,angle_min,0.0_MK)
           r2 = polar_star_r(a,b,c,angle_max,0.0_MK)
           
           x1 = r1*COS(angle_min)
           y1 = r1*SIN(angle_min)
           x2 = r2*COS(angle_max)
           y2 = r2*SIN(angle_max)
           
           d1 = SQRT((u-x1)**2+(v-y1)**2)
           d2 = SQRT((u-x2)**2+(v-y2)**2)
           
           IF ( d1<d2 ) THEN
              
              theta = angle_min
              d = d1
              r = r1
              
           ELSE
              
              theta = angle_max
              d = d2
              r = r2
              
           END IF
           
        ELSE IF ( ABS(f1) < mcf_machine_zero ) THEN
           
           r1 = polar_star_r(a,b,c,angle_min,0.0_MK)

           x1 = r1*COS(angle_min)
           y1 = r1*SIN(angle_min)

           d1 = SQRT((u-x1)**2+(v-y1)**2)
           
           theta = angle_min
           d     = d1
           r     = r1
           
           theta2 = &
                polar_star_F_root(a,b,c,u,v, &
                angle_min+mcf_machine_zero,&
                angle_max,stat_info_sub)
           
           IF ( theta2 >= 0.0_MK ) THEN
              
              r2 = polar_star_r(a,b,c,theta2,0.0_MK)
              
              x2 = r2*COS(theta2)
              y2 = r2*SIN(theta2)
              
              d2 = SQRT((u-x2)**2+(v-y2)**2)
              
              IF ( d2<d1 ) THEN
                 
                 theta = theta2
                 d = d2
                 r = r2
                 
              END IF
              
           END IF
           
        ELSE IF ( ABS(f2) < mcf_machine_zero ) THEN
           
           r2 = polar_star_r(a,b,c,angle_max,0.0_MK)
           
           x2 = r2*COS(angle_max)
           y2 = r2*SIN(angle_max)
           
           d2 = SQRT((u-x2)**2+(v-y2)**2)
           
           theta = angle_max
           d     = d2
           r     = r2
           
           theta1 = &
                polar_star_F_root(a,b,c,u,v, &
                angle_min,angle_max-mcf_machine_zero, &
                stat_info_sub)
           
           IF ( theta1 >= 0.0_MK ) THEN
              
              r1 = polar_star_r(a,b,c,theta1,0.0_MK)
              x1 = r1*COS(theta1)
              y1 = r1*SIN(theta1)
              d1 = SQRT((u-x1)**2+(v-y1)**2)
              
              IF ( d1<d2 ) THEN
                 
                 theta = theta1
                 d = d1
                 r = r1
                 
              END IF
              
           END IF
           
        ELSE
           
           theta = &
                polar_star_F_root(a,b,c,u,v, &
                angle_min,angle_max,stat_info_sub)
           
           r = polar_star_r(a,b,c,theta,0.0_MK)
           x = r*COS(theta)
           y = r*SIN(theta)
           d = SQRT((u-x)**2+(v-y)**2)
           
        END IF
        
        
        !----------------------------------------------------
        ! Flip the point from A(x,y) using symmetry.
        !
        ! Rotate counter-closewise (x,y) by phi radian.
        !----------------------------------------------------
        
        IF ( flip ) THEN
           
           theta = angle_sec - theta
           
        END IF
        
        theta = theta + angle_sec * sec
        
        theta = theta + phi
        
        x  = r*COS(theta)
        y  = r*SIN(theta)
        
        RETURN
        
      END SUBROUTINE polar_star_shortestD
      
      
    END MODULE COLLOID
    
    
      PROGRAM test
        
        USE COLLOID
        IMPLICIT NONE
        
        INTEGER                           :: num
        REAL(MK), DIMENSION(2)            :: sx
        REAL(MK), DIMENSION(2)            :: dx
        INTEGER                           :: id
        REAL(MK)                          :: a
        REAL(MK)                          :: b
        REAL(MK)                          :: c
        REAL(MK)                          :: d, d_short
        REAL(MK)                          :: r
        REAL(MK)                          :: theta
        REAL(MK), DIMENSION(2)            :: x


        num = 0
        
        a = 0.02_MK
        b = 0.01_MK
        c = 4.0_MK
        
        dx(:) = 2.0e-3_MK
        sx(2) = 0.0_MK + 0.5_MK * dx(2)
        
        
        DO WHILE( sx(2) <= 0.1_MK )
           
           sx(1) = 0.0_MK + 0.5_MK * dx(1)
           
           DO WHILE( sx(1) <= 0.1_MK )
              
              num = num + 1
              
              id = 0
              theta = angle(sx(1),sx(2))
              d = SQRT(sx(1)**2 + sx(2)**2)
              r = polar_star_r(a,b,c,theta,0.0_MK)
              
              IF ( d <= r) THEN
                 
                 id = 1
                 
                 CALL polar_star_shortestD(a,b,c,0.0_MK,sx(1),sx(2),&
                      x(1),x(2),d_short)
                 
                 !PRINT *, x(1),x(2), d_short
                 IF ( d_short > 6.0e-3 ) THEN
                    id = 2
                 END IF
                 
              END IF
              
              PRINT *, sx(1), sx(2), id
              
              sx(1) = sx(1) + dx(1)
              
           END DO ! sx(1)
           
           sx(2) = sx(2) + dx(2)
           
        END DO ! sx(2)
        
        
      END PROGRAM test
      
