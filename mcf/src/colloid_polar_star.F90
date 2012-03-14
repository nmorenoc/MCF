!----------------------------------------------------
! Star shapes related routines in polar coordinate system
!----------------------------------------------------
      REAL(MK) FUNCTION colloid_polar_star_r(a,b,c,theta,phi)
        !----------------------------------------------------
        ! Function value of 
        ! r(theta)=a+b*cos(c*theta-phi).
        !----------------------------------------------------

        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        colloid_polar_star_r = a + b * COS(c*(theta-phi))               
        
      END FUNCTION colloid_polar_star_r
      
      
      REAL(MK) FUNCTION colloid_polar_star_dr(b,c,theta,phi)
        !----------------------------------------------------
        ! Derivative value of 
        ! r(theta)=a+b*cos(c*theta-phi).
        !----------------------------------------------------

        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        colloid_polar_star_dr =  -b * c * SIN(c*(theta-phi))
        
      END FUNCTION colloid_polar_star_dr
      
      
      REAL(MK) FUNCTION colloid_polar_star_ddr(b,c,theta,phi)
        !----------------------------------------------------
        ! Second derivative value of 
        ! r(theta)=a+b*cos(c*theta-phi).
        !----------------------------------------------------

        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        colloid_polar_star_ddr =  -b * c * c * COS(c*(theta-phi))
        
      END FUNCTION colloid_polar_star_ddr
      

      REAL(MK) FUNCTION colloid_polar_star_G(a,b,c,t,u,v,p,q)
        !----------------------------------------------------
        ! Function G's value.
        ! Line connecting A(u,v) and B(p,q) is passing r(t)
        ! at C(x,y).
        ! Therefore, (y-q)(u-p)=(x-p)(v-q).
        ! Set G(t)=(v-q)x-(u-p)y+p(q-v)+q(u-p), where
        ! x=r(t)cos(t), y=r(t)sin(t) and r(t)=a+bcos(ct). 
        !----------------------------------------------------

        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: t
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        REAL(MK), INTENT(IN)            :: p
        REAL(MK), INTENT(IN)            :: q
        
        REAL(MK)                        :: r
        
        r = colloid_polar_star_r(a,b,c,t,0.0_MK)
        
        colloid_polar_star_G = (v-q)*r*COS(t) - &
             (u-p)*r*SIN(t) + p*(q-v)+q*(u-p)
        
        RETURN
        
      END FUNCTION colloid_polar_star_G
      
      
      REAL(MK) FUNCTION colloid_polar_star_dG(a,b,c,t,u,v,p,q)
        !----------------------------------------------------
        ! Function G's derivative value.
        ! Line connecting A(u,v) and B(p,q) is passing r(t)
        ! at C(x,y).
        ! Therefore, (y-q)(u-p)=(x-p)(v-q).
        ! Set G(t)=(v-q)x-(u-p)y+p(q-v)+q(u-p), where
        ! x=r(t)cos(t), y=r(t)sin(t) and r(t)=a+bcos(ct). 
        ! therefore, 
        ! G'(t) = (v-q)(r'cos(t)-rsin(t))
        !         -(u-p)(r'sin(t)+rcos(t))
        !----------------------------------------------------

        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: t
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        REAL(MK), INTENT(IN)            :: p
        REAL(MK), INTENT(IN)            :: q
        
        REAL(MK)                        :: r
        REAL(MK)                        :: dr
        
        
        r  = colloid_polar_star_r(a,b,c,t,0.0_MK)
        dr = colloid_polar_star_dr(b,c,t,0.0_MK)
        
        colloid_polar_star_dG = &
             (v-q)*(dr*COS(t) -r*SIN(t)) - &
             (u-p)*(dr*SIN(t)+ r*COS(t))
        
        RETURN
        
      END FUNCTION colloid_polar_star_dG
      
      
      SUBROUTINE zero_bracket_G(a,b,c,u,v,p,q,&
           x1,x2,n,xb1,xb2,nb,stat_info)
        !----------------------------------------------------
        ! Return
        ! nb: number of roots
        ! xb1,xb2:  bounds of roots.
        ! not including x1,x2.
        ! Function F's zero points lie inside [xb1,xb2].
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments:
        ! r(theta) = a * b*cos(c*theta).
        ! x1: theta_min
        ! x2: theta_max
        ! n : number of pieces
        !----------------------------------------------------

        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        REAL(MK), INTENT(IN)            :: p
        REAL(MK), INTENT(IN)            :: q     
        REAL(MK), INTENT(IN)            :: x1
        REAL(MK), INTENT(IN)            :: x2
        INTEGER, INTENT(IN)             :: n
        REAL(MK), DIMENSION(:), POINTER :: xb1
        REAL(MK), DIMENSION(:), POINTER :: xb2
        INTEGER, INTENT(INOUT)          :: nb
        INTEGER, INTENT(OUT)            :: stat_info

        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------

        INTEGER                         :: nn
        REAL(MK), DIMENSION(:), POINTER :: xt1
        REAL(MK), DIMENSION(:), POINTER :: xt2
        REAL(MK)                        :: x,gp,gc,dx
        INTEGER                         :: i,nbb
        
        !----------------------------------------------------
        ! Check input parameters and initialize variables.
        !----------------------------------------------------

        IF ( n < 1 ) THEN
           PRINT *, 'zero_bracket_G : ', &
                'interval n should be > 0 :', n
           stat_info = -1
           GOTO 9999
        END IF

        IF ( nb < 1 ) THEN
           PRINT *, 'zero_bracket_G : ', &
                'root nb should be > 0 :', nb
           stat_info = -1
           GOTO 9999
        END IF
        
        stat_info = 0
        
        NULLIFY(xt1)
        NULLIFY(xt2)
        
        ALLOCATE(xt1(n))
        ALLOCATE(xt2(n))
        
        !----------------------------------------------------
        ! Set resolution dx of searching roots.
        !----------------------------------------------------

        nbb = 0
        nn  = n
        IF ( ABS(x2-x1) <= mcf_machine_zero ) THEN
           nn = 1
        END IF
        
        dx = (x2-x1) / nn
        x  = x1
        gp = colloid_polar_star_G(a,b,c,x,u,v,p,q)
        
        DO i = 1, nn
           
           x = x + dx
           gc = colloid_polar_star_G(a,b,c,x,u,v,p,q)
           
           IF ( gc * gp <= mcf_machine_zero )  THEN
              
              nbb = nbb + 1
              xt1(nbb) = x -dx
              xt2(nbb) = x
              IF ( nbb == nb ) THEN 
                 GOTO 9999
              END IF
              
           END IF

           gp = gc
           
        END DO
        
        nb = nbb
        
9999    CONTINUE
        
        IF ( ASSOCIATED(xb1) ) THEN
           DEALLOCATE(xb1)
        END IF
        
        IF ( ASSOCIATED(xb2) ) THEN
           DEALLOCATE(xb2)
        END IF
        
        ALLOCATE(xb1(nb)) 
        ALLOCATE(xb2(nb))
        
        xb1(1:nb) = xt1(1:nb)
        xb2(1:nb) = xt2(1:nb)

        IF ( ASSOCIATED(xt1) ) THEN
           DEALLOCATE(xt1)
        END IF
        
        IF ( ASSOCIATED(xt2) ) THEN
           DEALLOCATE(xt2)
        END IF
        
        RETURN
        
      END SUBROUTINE zero_bracket_G
      
      
      REAL(MK) FUNCTION colloid_polar_star_G_root(a,b,c,u,v,p,q,&
           xb1,xb2,xacc,stat_info)
        
        !----------------------------------------------------
        ! Find the root of function G(t) = 0.0
        ! in (xb1, xb2).
        !
        ! Procedure :
        ! We use bisection+Newton-Raphson.
        ! See <<Numerical Recipes, Fortran version>>
        ! page 258.
        !----------------------------------------------------
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        REAL(MK), INTENT(IN)            :: p
        REAL(MK), INTENT(IN)            :: q       
        REAL(MK), INTENT(IN)            :: xb1
        REAL(MK), INTENT(IN)            :: xb2
        REAL(MK), INTENT(IN)            :: xacc
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here.
        !
        !----------------------------------------------------
        
        REAL(MK)                        :: x1,x2,xh,rts, root
        REAL(MK)                        :: dg,dx,dxold
        REAL(MK)                        :: g,gh,g1
        
        REAL(MK)                        :: temp
        
        INTEGER                         :: iter_max
        INTEGER                         :: iter
        
        !----------------------------------------------------
        ! Initialize variables.
        !----------------------------------------------------

        stat_info = 0
        iter_max = 100
        iter     = 0
        
        x1 = xb1
        x2 = xb2
        
        g1 = colloid_polar_star_G(a,b,c,x1,u,v,p,q)
        gh = colloid_polar_star_G(a,b,c,x2,u,v,p,q)
        
        IF ( g1*gh > ABS(mcf_machine_zero) ) THEN
           PRINT *, "colloid_polar_star_G_root : ", &
                "Root must be bracketed"
           PRINT *, "u, v : ", u,v
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( ABS(g1) <= mcf_machine_zero ) THEN
           root = x1
           GOTO 9999
        END IF
        
        IF ( ABS(gh) <= mcf_machine_zero ) THEN
           root = x2
           GOTO 9999
        END IF
        
        IF ( g1 < 0.0_MK ) THEN
           x1 = x1
           xh = x2
        ELSE
           xh = x1
           x1 = x2
        END IF
        
        rts   = (x1+x2) / 2.0_MK
        dxold = ABS(x2-x1)
        dx    = dxold
        
        g   = colloid_polar_star_G(a,b,c,rts,u,v,p,q)
        dg  = colloid_polar_star_dG(a,b,c,rts,u,v,p,q)
        
        DO iter = 1, iter_max
           
           IF ( ((rts-xh)*dg-g) * ( (rts-x1)*dg-g) > -mcf_machine_zero .OR. &
                ABS(2.0_MK*g) > ABS(dxold*dg)) THEN
              
              dxold = dx
              dx    = (xh-x1) / 2.0_MK
              rts   = x1+dx
              
              IF ( ABS(x1-rts) <= mcf_machine_zero ) THEN
                 
                 root = rts
                 GOTO 9999
                 
              END IF
              
           ELSE
              
              dxold = dx
              dx    = g/dg
              temp  = rts
              rts   = rts - dx
              
              IF ( ABS(temp-rts)<= mcf_machine_zero ) THEN
                 
                 root = rts
                 GOTO 9999
                 
              END IF
              
           END IF
           
           IF ( ABS(dx) <= xacc ) THEN
              
              root = rts 
              
           END IF
           
           IF ( g < 0.0_MK ) THEN
              
              x1 = rts
              
           ELSE
              
              xh = rts
              
           END IF
           
        END DO
        
        
9999    CONTINUE
        
        !PRINT *, "iter: ", iter
        
        IF ( iter >= iter_max ) THEN
           
           stat_info = -1
           
        ELSE
           
           colloid_polar_star_G_root = root
           
        END IF
        
        
        RETURN
        
      END FUNCTION colloid_polar_star_G_root
      
      
      SUBROUTINE colloid_polar_star_intersectP(a,b,c,phi,u,v,p,q,x,y,stat_info)
        !----------------------------------------------------
        ! Find the intersecting point of line A(u,v) and B(p,q)
        ! to the curve given in polar coordinate 
        ! r(t) = a+b*cos(c*(t-phi));
        ! Note that A(u,v) is inside curve, B(p,q) is outside,
        ! it is important.
        !
        ! Procedure :
        ! We first rotate axis by phi radian counter-clockwise,
        ! then we have standard star function a+bcos(ct).
        ! After calculationg, rotate back by phi radian.
        !
        ! Then adjust t_A, t_B to [0, 2pi).
        !
        ! 1: t_A<=t_B:
        !    a: t_B <= t_A+pi, search [t_A,t_B] for zero of G.
        !    b: t_B > t_A+pi, set t_B=t_B-2pi, 
        !       search [t_A, t_B] also for zero of G.
        ! 2: t_A >t_B::
        !    a: t_A < t_B+pi, search [t_A, t_B]
        !    b: t_A >= t_B+pi, set t_A=t_A-2pi,
        !       search [t_A, t_B].
        !----------------------------------------------------
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: phi
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v    
        REAL(MK), INTENT(IN)            :: p
        REAL(MK), INTENT(IN)            :: q
        REAL(MK), INTENT(OUT)           :: x
        REAL(MK), INTENT(OUT)           :: y
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here.
        !
        !----------------------------------------------------
        
        REAL(MK)                        :: thetaA,thetaB
        REAL(MK)                        :: tu,tv,tp,tq
        REAL(MK)                        :: rA,rB
        REAL(MK)                        :: theta,r
        REAL(MK), DIMENSION(:), POINTER :: xb1,xb2
        INTEGER                         :: nb       
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        NULLIFY(xb1)
        NULLIFY(xb1)
        nb = 1
        
        !----------------------------------------------------
        ! Calculate polar angles of A and B,
        ! check if A is inside and B is outside.
        !----------------------------------------------------
        
        thetaA = colloid_polar_angle(u,v)
        rA     = SQRT(u**2+v**2)
        r      = colloid_polar_star_r(a,b,c,thetaA, phi)
        
        IF ( rA > r ) THEN
           PRINT *, "colloid_polar_star_intersectP: ", &
                "A should be inside !"
           stat_info = -1
           GOTO 9999
        END IF
        
        thetaB =  colloid_polar_angle(p,q)
        rB     = SQRT(p**2+q**2)
        r      = colloid_polar_star_r(a,b,c,thetaB, phi)
        
        IF ( rB < r ) THEN
           PRINT *, "colloid_polar_star_intersectP: ", &
                "B should be outside !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! If (u,v) is at origin.
        !----------------------------------------------------

        IF ( ABS(u) <=mcf_machine_zero .AND. &
             ABS(v) <=mcf_machine_zero) THEN
           
           r  = colloid_polar_star_r(a,b,c,thetaB,0.0_MK)
           x  = r*COS(thetaB)
           y  = r*SIN(thetaB)
           GOTO 9999
           
        END IF
        
        !----------------------------------------------------
        ! Rotate axis by phi radian couter-clockwise.
        !----------------------------------------------------

        thetaA = thetaA - phi
        thetaB = thetaB - phi

        !----------------------------------------------------
        ! Adjust theta inbetween [0, 2pi).
        !----------------------------------------------------

        IF ( thetaA < 0.0_MK ) THEN
           thetaA = thetaA + 2.0_MK * mcf_pi
        ELSE IF ( thetaA >= 2.0_MK*mcf_pi) THEN
           thetaA = thetaA - 2.0_MK * mcf_pi
        END IF
        
        IF ( thetaB < 0.0_MK ) THEN
           thetaB = thetaB + 2.0_MK * mcf_pi
        ELSE IF ( thetaB >= 2.0_MK*mcf_pi) THEN
           thetaB = thetaB - 2.0_MK * mcf_pi
        END IF
        
        !----------------------------------------------------
        ! Check different situation, see the remarks in the
        ! head comment of this routine.
        !----------------------------------------------------
        
        !PRINT *, "thetaA,thetaB: ", thetaA, thetaB
        
        IF ( thetaA <= thetaB ) THEN
           
           IF ( thetaB > thetaA + mcf_pi ) THEN
              
              thetaB=thetaB-2.0_MK*mcf_pi
              
           END IF
           
        ELSE
           
           IF ( thetaA >= thetaB + mcf_pi ) THEN
              
              thetaA=thetaA-2.0_MK*mcf_pi
              
           END IF
           
        END IF
        
        !PRINT *, "thetaA,thetaB: ", thetaA, thetaB
        
        !----------------------------------------------------
        ! find point (tu,tv), (tp,tq) in new coordinate for
        ! (u,v) and (p,q).
        !----------------------------------------------------
        
        tu = rA*COS(thetaA)
        tv = rA*SIN(thetaA)
        
        tp = rB*COS(thetaB)
        tq = rB*SIN(thetaB)
        
        !PRINT *, thetaA, rA
        !PRINT *, thetaB, rB
         
        !PRINT *, "u,v,p,q: ", u,v,p,q
        !PRINT *, "tu,tv,tp,tq: ", tu,tv,tp,tq
        
        !----------------------------------------------------
        ! Find how many roots between two angles
        !  and bracket zeros into bounds.
        !----------------------------------------------------
        
        CALL zero_bracket_G(a,b,c,tu,tv,tp,tq,thetaB,thetaA,&
             50,xb1,xb2,nb,stat_info_sub)
        
        !PRINT *, "nb, xb1,xb2: ", nb, xb1(1:nb), xb2(1:nb)
        
        IF ( nb /= 1 ) THEN
           PRINT *, "colloid_polar_star_intersectP: nb != 1, wrong !"
           PRINT *, nb, xb1(nb), xb2(nb)
           stat_info = -1
           GOTO 9999
        END IF
        
        theta = &
             colloid_polar_star_G_root(a,b,c,tu,tv,tp,tq,xb1(1),xb2(1),&
             mcf_machine_zero*(xb1(1)+xb2(1))/2.0_MK,&
             stat_info_sub)
        
        IF ( stat_info_sub /= 0 )  THEN
           PRINT *, "colloid_polar_star_intersectP : ", &
                "No root"
           stat_info = -1
           GOTO 9999
        END IF
        
        r = colloid_polar_star_r(a,b,c,theta,0.0_MK)
        !PRINT *, theta, r
        !PRINT *
        
        !----------------------------------------------------
        ! Rotate counter-closewise (x,y) by phi radian.
        !----------------------------------------------------
        
        !print *, "phi: ", phi
        theta = theta + phi
        x  = r*COS(theta)
        y  = r*SIN(theta)
        
        !PRINT *, "x, y: ", x, y
        
        
9999    CONTINUE
        
        IF ( ASSOCIATED(xb1)) THEN
           DEALLOCATE(xb1)
        END IF
        
        IF ( ASSOCIATED(xb2)) THEN
           DEALLOCATE(xb2)
        END IF
        
        RETURN
        
      END SUBROUTINE colloid_polar_star_intersectP
      
      
      REAL(MK) FUNCTION colloid_polar_star_F(a,b,c,theta,u,v)
        !----------------------------------------------------
        ! Function F's value.
        !
        ! 1: tangent line's vector m of r(theta) is:
        !    a(dx/dt, dy/dt), i.e., 
        !    a(r'cos(t)-rsin(t), r'sin(t)+rcos(t))
        !    where a is a constant and t is angle theta.
        ! 2: vector n from (u,v) to tangent point is:
        !    (rcos(t)-u, rsin(t)-v)
        ! 3: obviously m (dot product) n = 0, therefore
        !    set F(t) = r*r'
        !               +r(usin(t)-vcos(t))
        !               -r'(vsin(t)+ucos(t)) = 0
        ! 4: Later on we solve equation F(t) = 0 to have
        !    t, then tangent point will be obtained.
        !----------------------------------------------------

        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        
        REAL(MK)                        :: r
        REAL(MK)                        :: dr
        
        !----------------------------------------------------
        ! Set last parameter phi to zero to caculate r and dr,
        ! since polar coordinate is assumed to be rotated phi
        ! counter-clockwise. 
        !----------------------------------------------------
        r  = colloid_polar_star_r(a,b,c,theta,0.0_MK)
        dr = colloid_polar_star_dr(b,c,theta,0.0_MK)
        
        colloid_polar_star_F = r*dr + &
             r  * ( u*SIN(theta) - v*COS(theta)) - &
             dr * ( v*SIN(theta) + u*COS(theta)) 
        
      END FUNCTION colloid_polar_star_F
      
      
      REAL(MK) FUNCTION colloid_polar_star_dF(a,b,c,theta,u,v)
        !----------------------------------------------------
        ! Above mentioned function F's derivative value.
        !----------------------------------------------------
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        
        REAL(MK)                        :: r
        REAL(MK)                        :: dr
        REAL(MK)                        :: ddr
        
        
        r  = colloid_polar_star_r(a,b,c,theta,0.0_MK)
        dr = colloid_polar_star_dr(b,c,theta,0.0_MK)
        ddr= colloid_polar_star_ddr(b,c,theta,0.0_MK)
        
        colloid_polar_star_dF = r*ddr + dr**2 + &
             2.0_MK*dr*( u*SIN(theta) - v*COS(theta)) + &
             ( r-ddr) *( u*COS(theta) + v*Sin(theta))
        
      END FUNCTION colloid_polar_star_dF

      
      SUBROUTINE zero_bracket_F(a,b,c,u,v, &
           x1,x2,n,xb1,xb2,nb,stat_info)
        !----------------------------------------------------
        ! Return
        ! nb: number of roots
        ! xb1,xb2:  bounds of roots.
        ! not including x1,x2.
        ! Function F's zero points lie inside [xb1,xb2].
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments:
        ! r(theta) = a * b*cos(c*theta).
        ! x1: theta_min
        ! x2: theta_max
        ! n : number of pieces
        !----------------------------------------------------

        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v        
        REAL(MK), INTENT(IN)            :: x1
        REAL(MK), INTENT(IN)            :: x2
        INTEGER, INTENT(IN)             :: n
        REAL(MK), DIMENSION(:), POINTER :: xb1
        REAL(MK), DIMENSION(:), POINTER :: xb2
        INTEGER, INTENT(INOUT)          :: nb
        INTEGER, INTENT(OUT)            :: stat_info

        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: nn
        REAL(MK), DIMENSION(:), POINTER :: xt1
        REAL(MK), DIMENSION(:), POINTER :: xt2
        REAL(MK)                        :: x,fp,fc,dx
        INTEGER                         :: i,nbb
        
        !----------------------------------------------------
        ! Check input parameters and initialize variables.
        !----------------------------------------------------

        IF ( n < 1 ) THEN
           PRINT *, 'zero_bracket_F : ', &
                'interval n should be > 0 :', n
           stat_info = -1
           GOTO 9999
        END IF

        IF ( nb < 1 ) THEN
           PRINT *, 'zero_bracket_F : ', &
                'root nb should be > 0 :', nb
           stat_info = -1
           GOTO 9999
        END IF
        
        stat_info = 0
        
        NULLIFY(xt1)
        NULLIFY(xt2)
        
        ALLOCATE(xt1(n))
        ALLOCATE(xt2(n))
        
        !----------------------------------------------------
        ! Set resolution dx of searching roots.
        !----------------------------------------------------

        nbb = 0
        nn  = n
        IF ( ABS(x2-x1)<=mcf_machine_zero ) THEN
           nn = 1
        END IF
        
        dx = (x2-x1) / nn
        x  = x1
        fp = colloid_polar_star_F(a,b,c,x,u,v)
        
        DO i = 1, nn

           x  = x + dx
           fc = colloid_polar_star_F(a,b,c,x,u,v)
           
           IF ( fc * fp <= mcf_machine_zero )  THEN
              
              nbb = nbb + 1
              xt1(nbb) = x -dx
              xt2(nbb) = x
              IF ( nbb == nb ) THEN 
                 GOTO 9999
              END IF
              
           END IF

           fp = fc
           
        END DO

        nb = nbb
        
9999    CONTINUE

        IF ( ASSOCIATED(xb1) ) THEN
           DEALLOCATE(xb1)
        END IF
        
        IF ( ASSOCIATED(xb2) ) THEN
           DEALLOCATE(xb2)
        END IF
        
        ALLOCATE(xb1(nb)) 
        ALLOCATE(xb2(nb))
        
        xb1(1:nb) = xt1(1:nb)
        xb2(1:nb) = xt2(1:nb)

        IF ( ASSOCIATED(xt1) ) THEN
           DEALLOCATE(xt1)
        END IF
        
        IF ( ASSOCIATED(xt2) ) THEN
           DEALLOCATE(xt2)
        END IF
        
        RETURN
        
      END SUBROUTINE zero_bracket_F
      
      
      REAL(MK) FUNCTION colloid_polar_star_F_root(a,b,c,u,v,&
           xb1,xb2,xacc,stat_info)
        
        !----------------------------------------------------
        ! Find the root of function F(theta) = 0.0
        ! in (theta_min, theta_max).
        !
        ! Procedure :
        ! We use bisection+Newton-Raphson.
        ! See <<Numerical Recipes, Fortran version>>
        ! page 258.
        !----------------------------------------------------
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: u
        REAL(MK), INTENT(IN)            :: v
        REAL(MK), INTENT(IN)            :: xb1
        REAL(MK), INTENT(IN)            :: xb2
        REAL(MK), INTENT(IN)            :: xacc
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here.
        !
        !----------------------------------------------------
   
        REAL(MK)                        :: x1,x2,xh,rts, root
        REAL(MK)                        :: df,dx,dxold
        REAL(MK)                        :: f,fh,f1
        
        REAL(MK)                        :: temp
        
        INTEGER                         :: iter_max
        INTEGER                         :: iter
        
        !----------------------------------------------------
        ! Initialize variables.
        !----------------------------------------------------

        stat_info = 0
        iter_max = 100
        iter     = 0

        x1 = xb1
        x2 = xb2
        
        f1 = colloid_polar_star_F(a,b,c,x1,u,v)
        fh = colloid_polar_star_F(a,b,c,x2,u,v)
        
        IF ( f1*fh > ABS(mcf_machine_zero) ) THEN
           PRINT *, "colloid_polar_star_F_root : ", &
                "Root must be bracketed"
           PRINT *, "u, v : ", u,v
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( ABS(f1) <= mcf_machine_zero ) THEN
           root = x1
           GOTO 9999
        END IF
        
        IF ( ABS(fh) <= mcf_machine_zero ) THEN
           root = x2
           GOTO 9999
        END IF
        
        IF ( f1 < 0.0_MK ) THEN
           x1 = x1
           xh = x2
        ELSE
           xh = x1
           x1 = x2
        END IF
        
        rts   = (x1+x2) / 2.0_MK
        dxold = ABS(x2-x1)
        dx    = dxold
        
        f   = colloid_polar_star_F(a,b,c,rts,u,v)
        df  = colloid_polar_star_dF(a,b,c,rts,u,v)
        
        DO iter = 1, iter_max
           
           IF ( ((rts-xh)*df-f) * ( (rts-x1)*df-f) > -mcf_machine_zero .OR. &
                ABS(2.0_MK*f) > ABS(dxold*df)) THEN
              
              dxold = dx
              dx    = (xh-x1) / 2.0_MK
              rts   = x1+dx
              
              IF ( ABS(x1-rts) <= mcf_machine_zero ) THEN
                 
                 root = rts
                 GOTO 9999
                 
              END IF
              
           ELSE
              
              dxold = dx
              dx    = f/df
              temp  = rts
              rts   = rts - dx
              
              IF ( ABS(temp-rts)<= mcf_machine_zero ) THEN
                 
                 root = rts
                 GOTO 9999
                 
              END IF
              
           END IF
           
           IF ( ABS(dx) < xacc ) THEN
              
              root = rts 
              
           END IF
           
           IF ( f < 0.0_MK ) THEN
              
              x1 = rts
              
           ELSE
              
              xh = rts
              
           END IF
           
        END DO
        
        
9999    CONTINUE
        
        ! PRINT *, iter
        IF ( iter >= iter_max ) THEN
           
           colloid_polar_star_F_root = -1.0_MK
           
        ELSE
           
           colloid_polar_star_F_root = root
           
        END IF
        
        
        RETURN
        
      END FUNCTION colloid_polar_star_F_root
      
      
      SUBROUTINE colloid_polar_star_shortestD(a,b,c,phi,p,q,x,y,d,stat_info)
        !----------------------------------------------------
        ! Find the shortest distance between a point (p,q)
        ! to the curve given in polar coordinate 
        ! r(t) = a+b*cos(c*(t-phi));
        ! Due to symmetry, single root is always found,
        ! i.e., no double/triple... roots.
        !
        ! Procedure :
        ! We first rotate axis by phi radian counter-clockwise,
        ! then always map point(p,q) by symmetry axis to
        ! the t=[0:pi/c] and do calculation there
        ! to find x,y,d by using polar_star_F_root
        ! (bisection).
        !
        ! Afterwards map (x,y) back and 
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
        INTEGER, INTENT(OUT)            :: stat_info
        
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
        REAL(MK)                        :: r1
        REAL(MK)                        :: x1,y1
        REAL(MK)                        :: d1
        REAL(MK), DIMENSION(:), POINTER :: xb1,xb2
        INTEGER                         :: nb, i
        
        REAL(MK)                        :: theta1

        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        d = 0.0_MK
        d = HUGE(d)
        flip = .FALSE.
        
        NULLIFY(xb1)
        NULLIFY(xb1)
        nb = 3
       
        !----------------------------------------------------
        ! Rotate axis by phi radian couter-clockwise.
        !----------------------------------------------------
        
        theta = colloid_polar_angle(p,q)
        r = SQRT(p**2+q**2)        
        theta = theta - phi
        
        !----------------------------------------------------
        ! Adjust theta inbetween [0, 2pi).
        !----------------------------------------------------

        IF ( theta < 0.0_MK ) THEN
           theta = theta + 2.0_MK * mcf_pi
        ELSE IF ( theta >= 2.0_MK*mcf_pi) THEN
           theta = theta - 2.0_MK * mcf_pi
        END IF
        
        !----------------------------------------------------
        ! Rotate theta to first section of rotation symmetry.
        !----------------------------------------------------
        
        angle_sec = 2.0_MK * mcf_pi / c
        sec       = INT(theta/angle_sec)
        theta     = theta - angle_sec * sec
        
        !----------------------------------------------------
        ! Flip angle in first section due to mirror symmetry.
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
        ! Find how many roots between angle_min and angle_max
        ! and bracket them into bounds.
        !----------------------------------------------------
        
        CALL zero_bracket_F(a,b,c,u,v,angle_min,angle_max,&
             50,xb1,xb2,nb,stat_info_sub)
        
        !PRINT *, nb, xb1(nb), xb2(nb)
        
        DO i = 1, nb
           
           theta1 = &
                colloid_polar_star_F_root(a,b,c,u,v,xb1(i),xb2(i),&
                mcf_machine_zero*(xb1(i)+xb2(i))/2.0_MK,&
                stat_info_sub)
           
           IF ( theta1 < 0.0_MK )  THEN
              PRINT *, "polar_star_shortestD : ", &
                   "No root"
              stat_info = -1
              GOTO 9999
           END IF
           
           r1 = colloid_polar_star_r(a,b,c,theta1,0.0_MK)
           x1 = r1*COS(theta1)
           y1 = r1*SIN(theta1)
           d1 = SQRT((u-x1)**2+(v-y1)**2)
           
           IF ( d1 < d ) THEN
              
              theta = theta1
              d     = d1
              r     = r1
              
           END IF
           
        END DO
        
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
        
        
9999    CONTINUE
        
        IF ( ASSOCIATED(xb1)) THEN
           DEALLOCATE(xb1)
        END IF
        
        IF ( ASSOCIATED(xb2)) THEN
           DEALLOCATE(xb2)
        END IF
        
        RETURN
        
      END SUBROUTINE colloid_polar_star_shortestD
      
      
      LOGICAL FUNCTION colloid_polar_star_convex(c,phi,p,q,stat_info)
        !----------------------------------------------------
        ! We first rotate axis by phi radian counter-clockwise,
        ! then always map point(p,q) by symmetry axis to
        ! the t=[0:pi/c].
        ! Then we check if it lies in the convex or concave
        ! part of the curve.
        ! The logic of this routine is only suitable for
        ! r(t)=a+bcos(ct) shape stars.
        !----------------------------------------------------
        
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: phi
        REAL(MK), INTENT(IN)            :: p
        REAL(MK), INTENT(IN)            :: q
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here.
        !
        ! angle_sec : angle of each section of symmetry.
        !----------------------------------------------------
        
        REAL(MK)                        :: theta
        REAL(MK)                        :: r
        REAL(MK)                        :: angle_sec
        INTEGER                         :: sec
        
        stat_info = 0
        
        !----------------------------------------------------
        ! Rotate axis by phi radian couter-clockwise.
        !----------------------------------------------------
        
        theta = colloid_polar_angle(p,q)
        r = SQRT(p**2+q**2)
        theta = theta - phi
        
        !----------------------------------------------------
        ! Adjust theta inbetween [0, 2pi).
        !----------------------------------------------------
        
        IF ( theta < 0.0_MK ) THEN
           theta = theta + 2.0_MK * mcf_pi
        ELSE IF ( theta >= 2.0_MK*mcf_pi ) THEN
           theta = theta - 2.0_MK * mcf_pi
        END IF
        
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
           
           theta = angle_sec - theta
           
        END IF
        
        IF ( theta <= mcf_pi / c / 2.0_MK ) THEN
           
           colloid_polar_star_convex = .TRUE.
           
        ELSE
           
           colloid_polar_star_convex = .FALSE.
           
        END IF
        
        RETURN
        
      END FUNCTION colloid_polar_star_convex
