      MODULE header
        
        IMPLICIT NONE

        !----------------------------------------------------
        ! Commonly used data.
        !----------------------------------------------------
        
        INTEGER, PARAMETER                      :: MK= KIND(1.0D0)
        
        REAL(MK),DIMENSION(:,:),POINTER         :: comb_value
        LOGICAL,DIMENSION(:,:), POINTER         :: comb_mark
        
        REAL(MK),DIMENSION(:,:,:),POINTER       :: p_value
        LOGICAL,DIMENSION(:,:,:),POINTER        :: p_mark
        REAL(MK),DIMENSION(:,:,:),POINTER       :: v_value
        LOGICAL,DIMENSION(:,:,:),POINTER        :: v_mark
        REAL(MK),DIMENSION(:,:,:),POINTER       :: q_value
        LOGICAL,DIMENSION(:,:,:),POINTER        :: q_mark
        
      CONTAINS
        
        
        REAL(MK) FUNCTION combination(n,k)
          !--------------------------------------------------
          ! Calcualting binominals/combinations.
          !--------------------------------------------------
          
          INTEGER, INTENT(IN)           :: n
          INTEGER, INTENT(IN)           :: k
          
          REAL(MK)                      :: c          
          INTEGER                       :: m
          INTEGER                       :: i
          
          IF ( k > n ) THEN
             PRINT *, "combination : k > n !"
             STOP
          END IF
          
          IF ( n-k < k ) THEN
             m = n-k
          ELSE
             m = k
          END IF
          
          c = 1
          
          DO i = 0, m-1
             c = c *(n-i)
          END DO
          
          DO i = 1, m
             c = c / i
          END DO
          
          combination = c
          
          RETURN
          
        END FUNCTION combination
        
        
        RECURSIVE REAL(MK) FUNCTION p_func(n,p,q)
          
          INTEGER, INTENT(IN)           :: n
          INTEGER, INTENT(IN)           :: p
          INTEGER, INTENT(IN)           :: q
          
          REAL(MK)                      :: pp
          REAL(MK)                      :: tp
          REAL(MK)                      :: t0,t1,t2,t3,t4
          INTEGER                       :: s
          
          pp = 0.0_MK
          
          IF ( p == 0 .AND. q == 0 ) THEN
             
             IF ( n == 1 ) THEN
                
                pp = 1.0_MK

             ELSE
                
                pp = 0.0_MK
                
             END IF
             
          ELSE
             
             t0 = 2*n+1
             t0 = t0/2/(n+1)
             t2 = n*(2*n-1)
             t2 = t2/2/(n+1)
             t4 = 2*(4*n**2-1)
             t4 = t4/3/(n+1)
             
             DO s = 1, q
                
                t1 = 3*(n+s)-(n*s+1)*(2*n*s-s-n+2)
                t1 = t1/s/(n+s)/(2*s-1)
                t1 = t0*t1
                t3 = n*(4*n**2-1)
                t3 = t3/2/(n+1)/(2*s+1)

                tp = 0.0_MK
                
                IF ( p_mark(s,q-s,p-n+1) ) THEN
                   tp = tp + t1 * p_value(s,q-s,p-n+1)
                ELSE
                   tp = tp + t1 * p_func(s,q-s,p-n+1)
                END IF
                
                IF ( p_mark(s,q-s,p-n-1) ) THEN
                   tp = tp + t2 * p_value(s,q-s,p-n-1)
                ELSE
                   tp = tp + t2 * p_func(s,q-s,p-n-1)
                END IF
                
                IF ( v_mark(s,q-s-2,p-n+1) ) THEN
                   tp = tp + t3 * v_value(s,q-s-2,p-n+1)
                ELSE
                   tp = tp + t3 * v_func(s,q-s-2,p-n+1)
                END IF
                
                IF ( q_mark(s,q-s-1,p-n+1) ) THEN
                   tp = tp - t4 * q_value(s,q-s-1,p-n+1)
                ELSE
                   tp = tp - t4 * q_func(s,q-s-1,p-n+1)
                END IF

                IF ( comb_mark(n+s,n+1) ) THEN
                   pp = pp + tp * comb_value(n+s,n+1)
                ELSE
                   pp = pp + tp * combination(n+s,n+1)
                END IF
                
             END DO
             
          END IF
          
          p_func = pp
          
          !----------------------------------------
          ! Record p_value globally and mark
          ! it as recorded.
          !----------------------------------------
          
          p_value(n,p,q) = pp
          p_mark(n,p,q)  = .TRUE.
          
          RETURN
          
        END FUNCTION p_func
        
        
        RECURSIVE REAL(MK) FUNCTION v_func(n,p,q)
          
          INTEGER, INTENT(IN)           :: n
          INTEGER, INTENT(IN)           :: p
          INTEGER, INTENT(IN)           :: q
          
          REAL(MK)                      :: vv
          REAL(MK)                      :: tv
          INTEGER                       :: s
          
          
          vv = 0.0_MK
          
          IF ( p == 0 .AND. q == 0 ) THEN
             
             IF ( n==1 ) THEN
                
                vv = 1.0_MK
                
             ELSE
                
                vv = 0.0_MK
                
             END IF
             
          ELSE
             
             DO s = 1,q
                
                tv = 0.0_MK
                
                IF ( p_mark(s,q-s,p-n+1) ) THEN
                   tv =  p_value(s,q-s,p-n+1)
                ELSE
                   tv =  p_func(s,q-s,p-n+1)
                END IF
                
                IF ( comb_mark(n+s,n+1) ) THEN
                   tv = tv * comb_value(n+s,n+1)
                ELSE
                   tv = tv * combination(n+s,n+1)
                END IF
                
                vv = vv + tv
                
             END DO
             
             vv = vv * 2*n /(n+1)/(2*n+3)
             
             IF ( p_mark(n,p,q) ) THEN
                vv = vv + p_value(n,p,q)
             ELSE
                vv = vv + p_func(n,p,q)
             END IF
             
          END IF
          
          v_func = vv
          
          !----------------------------------------
          ! Record v_value globally and mark
          ! it as recorded.
          !----------------------------------------
          
          v_value(n,p,q) = vv
          v_mark(n,p,q)  = .TRUE.
          
        END FUNCTION v_func
        

        RECURSIVE REAL(MK) FUNCTION q_func(n,p,q)
          
          INTEGER, INTENT(IN)           :: n
          INTEGER, INTENT(IN)           :: p
          INTEGER, INTENT(IN)           :: q
          
          REAL(MK)                      :: qq
          REAL(MK)                      :: tq
          INTEGER                       :: s
          
          
          qq = 0.0_MK
          
          IF ( p == 0 .AND. q == 0 ) THEN
             
             qq = 0.0_MK
             
          ELSE
             
             DO s = 1,q
             
                tq = 0.0_MK
                
                IF ( q_mark(s,q-s-1,p-n) ) THEN
                   tq = tq + &
                        q_value(s,q-s-1,p-n)*s/(n+1)
                ELSE
                   tq = tq + &
                        q_func(s,q-s-1,p-n)*s/(n+1)
                END IF
                
                IF ( p_mark(s,q-s,p-n) ) THEN
                   tq = tq - &
                        p_value(s,q-s,p-n)*3/2/n/s/(n+1)
                ELSE
                   tq = tq - &
                        p_func(s,q-s,p-n)*3/2/n/s/(n+1)
                END IF
                
                IF ( comb_mark(n+s,n+1) ) THEN
                   qq = qq + comb_value(n+s,n+1) * tq
                ELSE
                   qq = qq + combination(n+s,n+1) * tq
                END IF
                
             END DO
             
          END IF
          
          q_func = qq

          !----------------------------------------
          ! Record v_value globally and mark
          ! it as recorded.
          !----------------------------------------
          
          q_value(n,p,q) = qq
          q_mark(n,p,q)  = .TRUE.
          
        END FUNCTION q_func
        
        
      END MODULE header
      
      
      PROGRAM resistance
        
        !----------------------------------------------------
        ! Main routine :
        !
        ! Resistance problem of shear movement
        ! complete/full calculation
        ! according to Jeffrey and Onishi 1984.
        !----------------------------------------------------
        
        USE header
        
        IMPLICIT NONE
        
        !----------------------------------------------------
        ! a1,a2,a : radius of two spheres
        ! R       : initial distance between centers
        ! b       : a2/a1
        !----------------------------------------------------
        
        REAL(MK),PARAMETER              :: pi=3.14159265358979_MK
        REAL(MK)                        :: a
        REAL(MK)                        :: R
        REAL(MK)                        :: d,dd
        REAL(MK)                        :: lambda
        REAL(MK),DIMENSION(:),POINTER   :: f
        INTEGER                         :: k
        REAL(MK)                        :: eta
        REAL(MK)                        :: YA11
        REAL(MK)                        :: YA12
        REAL(MK)                        :: U
        REAL(MK)                        :: Fd
        CHARACTER(256)                  :: cbuf
        INTEGER                         :: n_arg
        INTEGER                         :: n,p,q,s
        REAL(MK)                        :: temp
        
        REAL(MK)                        :: s_min, s_max
        INTEGER                         :: num_p
        REAL(MK)                        :: s_cur,s_dif
        INTEGER                         :: num_k
        
        !------------------------------------------
        ! Initialization of variables.
        !------------------------------------------
        
        R = 0.08_MK
        a = 0.02_MK
        lambda = 1.0_MK
        
        s_min = 2.0_MK
        s_max = 10.0_MK
        num_p = 80
        num_k = 12
        
        eta  = 0.1_MK
        YA11 = 0.0_MK
        YA12 = 0.0_MK
        
        U    = 1.25e-4_MK
        Fd   = 0.0_MK
        
        NULLIFY(comb_value)
        NULLIFY(comb_mark)
        NULLIFY(p_value)
        NULLIFY(p_mark)
        NULLIFY(v_value)
        NULLIFY(v_mark)
        NULLIFY(q_value)
        NULLIFY(q_mark)
        
        NULLIFY(f)
        
        n_arg = IARGC()
        
        !------------------------------------------
        ! s_min : minimal gap.
        ! s_max : maximal gap.
        ! num_p : number of points in gap.
        ! num_k : number of series points for accuracy.
        !------------------------------------------
        
        IF ( n_arg > 3 ) THEN
           CALL GETARG(1,cbuf)
           READ(cbuf,*), s_min
           CALL GETARG(2,cbuf)
           READ(cbuf,*), s_max
           CALL GETARG(3,cbuf)
           READ(cbuf,*), num_p         
           CALL GETARG(4,cbuf)
           READ(cbuf,*), num_k
        END IF
        
        !------------------------------------------
        ! number of series points are always even
        ! number started from 0 -> num_k
        !------------------------------------------
        
        IF ( MOD(num_k,2) /= 0 ) THEN
           num_k = num_k + 1
        END IF
        
        !------------------------------------------
        ! distance between points in the gap.
        !------------------------------------------
        
        s_dif = (s_max-s_min)/num_p
        
        !------------------------------------------
        ! Allocate memory to record combination.
        !------------------------------------------
        
        ALLOCATE(comb_value(0:2*num_k,0:num_k+1))
        ALLOCATE(comb_mark(0:2*num_k,0:num_k+1))

        comb_value(:,:) = 0.0_MK
        comb_mark(:,:)  = .FALSE.

        !------------------------------------------
        ! Calculate all combinations/binomials needed.
        !------------------------------------------

        DO n = -1, num_k
           DO s = 1, num_k
              comb_value(n+s,n+1) = combination(n+s,n+1)
              comb_mark(n+s,n+1) = .TRUE.
           END DO
        END DO
        
        !------------------------------------------
        ! Allocate memory to record p and v values,
        ! and their marks.e
        !------------------------------------------
        
        ALLOCATE(p_value(-num_k-1:num_k,-num_k-1:num_k,-num_k-1:num_k))
        ALLOCATE(p_mark (-num_k-1:num_k,-num_k-1:num_k,-num_k-1:num_k))
        ALLOCATE(v_value(-num_k-1:num_k,-num_k-1:num_k,-num_k-1:num_k))
        ALLOCATE(v_mark (-num_k-1:num_k,-num_k-1:num_k,-num_k-1:num_k))
        ALLOCATE(q_value(-num_k-1:num_k,-num_k-1:num_k,-num_k-1:num_k))
        ALLOCATE(q_mark (-num_k-1:num_k,-num_k-1:num_k,-num_k-1:num_k))
     
        p_value(:,:,:) = 0.0_MK
        p_mark(:,:,:)  = .FALSE.
        v_value(:,:,:) = 0.0_MK
        v_mark(:,:,:)  = .FALSE.
        q_value(:,:,:) = 0.0_MK
        q_mark(:,:,:)  = .FALSE.
        
        !------------------------------------------
        ! Allocate memory for f function values.
        !------------------------------------------
        
        ALLOCATE(f(0:num_k-1))
        f(:) = 0.0_MK
        
        !------------------------------------------
        ! coefficients.
        !------------------------------------------
        
        DO k = 0, num_k-1
           
           DO q = 0, k
              
              IF ( p_mark(1,k-1,q) ) THEN
                 f(k) = f(k) + p_value(1,k-q,q) * lambda**q
                 PRINT *, 2**k*p_value(1,k-q,q)                 
              ELSE
                 f(k) = f(k) + p_func(1,k-q,q)  * lambda**q
                 PRINT *, 2**k*p_func(1,k-q,q)
              END IF
              
           END DO
           
           f(k) = f(k) * 2**k
           !PRINT *, f(k)
           PRINT *, ""
           
        END DO
        
#if 0        
        f(0) = 1.0_MK
        f(1) = 1.5_MK
        f(2) = 9.0_MK/4.0_MK
        f(3) = 4.0_MK + 27.0_MK / 8.0_MK 
        f(4) = 24.0_MK + 81.0_MK /16.0_MK
        f(5) = 63.0_MK + 243.0_MK/32.0_MK
        f(6) = 211.0_MK + 1241.0_MK/64.0_MK
        f(7) = 288.0_MK + 1053.0_MK/8.0_MK + &
             19083.0_MK/128.0_MK + 1053.0_MK/8.0_MK
        f(8) = 1215.0_MK + 4261.0_MK/8.0_MK + &
             126369.0_MK/256.0_MK - 117.0_MK/8.0_MK 
        f(9) = 1152.0_MK + 2268.0_MK + &
             60443.0_MK /16.0_MK + 766179.0_MK/512.0_MK
        f(10) = 2304.0_MK + 7857.0_MK/4.0_MK + &
             98487.0_MK/16.0_MK + 10548393.0_MK/1024.0_MK + &
             67617.0_MK/8.0_MK - 351.0_MK/2.0_MK + 3888.0_MK
        f(11) = 4608.0_MK + 14256.0_MK + 22071.0_MK + &
             2744505.0_MK/64.0_MK + 95203835.0_MK/2048.0_MK
#endif
        
        DO p=0, num_p
           
           s_cur = s_min + p * s_dif
           
           d = 1.0_MK/s_cur/(1.0_MK+lambda)
           
           dd = 1.0_MK
           
           YA11 = 0.0_MK
           YA12 = 0.0_MK
           
           !PRINT *, d, dd, (num_k-1)/2
           DO k = 0, (num_k-1)/2
              
              YA11 = YA11 + f(2*k) * dd
              
              YA12 = YA12 - f(2*k+1) * dd * d
              
              dd = dd*d**2
              
           END DO
           
           
           Fd = 6.0_MK*pi*a*eta*(YA11-YA12)*U
           !Fd = (YA11-YA12)
           
           PRINT *, s_cur-2.0_MK, Fd
           
        END DO
        
9999    CONTINUE
        
        IF ( ASSOCIATED(comb_value)) THEN
           DEALLOCATE(comb_value)
        END IF

        IF ( ASSOCIATED(comb_mark)) THEN
           DEALLOCATE(comb_mark)
        END IF
        
        IF ( ASSOCIATED(p_value)) THEN
           DEALLOCATE(p_value)
        END IF
        
        IF ( ASSOCIATED(p_mark)) THEN
           DEALLOCATE(p_mark)
        END IF
        
        IF ( ASSOCIATED(v_value)) THEN
           DEALLOCATE(v_value)
        END IF
        
        IF ( ASSOCIATED(v_mark)) THEN
           DEALLOCATE(v_mark)
        END IF

        IF ( ASSOCIATED(q_value)) THEN
           DEALLOCATE(q_value)
        END IF
        
        IF ( ASSOCIATED(q_mark)) THEN
           DEALLOCATE(q_mark)
        END IF
        
        IF ( ASSOCIATED(f)) THEN
           DEALLOCATE(f)
        END IF
        
        
      END PROGRAM resistance
       
      
