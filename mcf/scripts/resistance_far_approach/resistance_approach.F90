      MODULE header
        
        IMPLICIT NONE

        !------------------------------------------
        ! Commonly used data.
        !------------------------------------------
        
        INTEGER, PARAMETER                      :: MK= KIND(1.0D0)
        
        REAL(MK),DIMENSION(:,:),POINTER         :: comb_value
        
        REAL(MK),DIMENSION(:,:,:),POINTER       :: p_value
        LOGICAL,DIMENSION(:,:,:),POINTER        :: p_mark
        REAL(MK),DIMENSION(:,:,:),POINTER       :: v_value
        LOGICAL,DIMENSION(:,:,:),POINTER        :: v_mark
        
      CONTAINS
        
        REAL(MK) FUNCTION combination(n,k)
          
          INTEGER, INTENT(IN)           :: n
          INTEGER, INTENT(IN)           :: k

          INTEGER                       :: m
          REAL(MK)                      :: c
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
          REAL(MK)                      :: t1,t2,t3 
          INTEGER                       :: s
          
          pp = 0.0_MK
          
          IF ( p == 0 .AND. q == 0 ) THEN
             
             IF ( n == 1 ) THEN
                
                pp = 1.0_MK

             ELSE
                
                pp = 0.0_MK
                
             END IF
             
          ELSE
             
             t2 = n*(2*n-1)
             t2 = t2/2/(n+1)
             pp = 0.0_MK
             
             DO s = 1, q
              
                t1 = n*(2*n+1)*(2*n*s-n-s+2)
                t1 = t1/2/(n+1)/(2*s-1)/(n+s)
                t3 = REAL(n,MK)*(4*n**2-1)/2/(n+1)/(2*s+1)

                tp = 0.0_MK
                
                !tp = t1 * p_func(s,q-s,p-n+1) - &
                !     t2 * p_func(s,q-s,p-n-1) - &
                !     t3 * v_func(s,q-s-2,p-n+1)
                
                IF ( p_mark(s,q-s,p-n+1) ) THEN
                   tp = tp + t1 * p_value(s,q-s,p-n+1)
                ELSE
                   tp = tp + t1 * p_func(s,q-s,p-n+1)
                END IF
                
                IF ( p_mark(s,q-s,p-n-1) ) THEN
                   tp = tp - t2 * p_value(s,q-s,p-n-1)
                ELSE
                   tp = tp - t2 * p_func(s,q-s,p-n-1)
                END IF
                
                IF ( v_mark(s,q-s-2,p-n+1) ) THEN
                   tp = tp - t3 * v_value(s,q-s-2,p-n+1)
                ELSE
                   tp = tp - t3 * v_func(s,q-s-2,p-n+1)
                END IF
                
                !pp = pp + tp * combination(n+s,n)
                pp = pp + tp * comb_value(n+s,n)
                
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
             
             tv = 0.0_MK
             
             DO s = 1,q
                
                !PRINT *, "P(", s, q-s-s,p-n-1, ")"
                !PRINT *, ""
                
                !tv = tv + &
                !     p_func(s,q-s,p-n-1) * combination(n+s,n)
                
                IF ( p_mark(s,q-s,p-n-1) ) THEN
                   tv = tv + &
                        p_value(s,q-s,p-n-1) * comb_value(n+s,n)
                   
                ELSE
                   tv = tv + &
                        p_func(s,q-s,p-n-1) * combination(n+s,n)
                END IF
                
             END DO
             
             vv = p_func(n,p,q) - &
                  tv*2*n /(n+1)/(2*n+3)
             
             !PRINT *, "P(", n, p,q, ")"
             !PRINT *, ""
             
          END IF
          
          v_func = vv

          !----------------------------------------
          ! Record v_value globally and mark
          ! it as recorded.
          !----------------------------------------
          
          v_value(n,p,q) = vv
          v_mark(n,p,q)  = .TRUE.
          
        END FUNCTION v_func
        
      END MODULE header
      
      
      PROGRAM resistance
        
        !----------------------------------------------------
        ! Main routine :
        !
        ! Resistance problem complete/full calculation
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
        REAL(MK),DIMENSION(0:10)        :: b
        REAL(MK)                        :: R
        REAL(MK)                        :: d,dd
        REAL(MK)                        :: lambda
        REAL(MK),DIMENSION(:),POINTER   :: f
        INTEGER                         :: k
        REAL(MK)                        :: eta
        REAL(MK)                        :: XA11
        REAL(MK)                        :: XA12
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
        s_max = 5.0_MK
        num_p = 30
        num_k = 12
        
        eta  = 0.1_MK
        XA11 = 0.0_MK
        XA12 = 0.0_MK
        
        U    = 1.25e-4_MK
        Fd   = 0.0_MK
        
        NULLIFY(comb_value)        
        NULLIFY(p_value)
        NULLIFY(p_mark)
        NULLIFY(v_value)
        NULLIFY(v_mark)
        
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
        
        ALLOCATE(comb_value(0:2*num_k-2,0:num_k-1))
        
        comb_value(:,:) = 0.0_MK
        
        !------------------------------------------
        ! Allocate memory to record p and v values,
        ! and their marks.e
        !------------------------------------------
        
        ALLOCATE(p_value(-num_k-1:num_k,-num_k-1:num_k,-num_k-1:num_k))
        ALLOCATE(p_mark(-num_k-1:num_k,-num_k-1:num_k,-num_k-1:num_k))
        ALLOCATE(v_value(-num_k-1:num_k,-num_k-1:num_k,-num_k-1:num_k))
        ALLOCATE(v_mark(-num_k-1:num_k,-num_k-1:num_k,-num_k-1:num_k))
        
        p_value(:,:,:) = 0.0_MK
        p_mark(:,:,:)  = .FALSE.
        v_value(:,:,:) = 0.0_MK
        v_mark(:,:,:)  = .FALSE.
        
        !------------------------------------------
        ! Allocate memory for f function values.
        !------------------------------------------
        
        ALLOCATE(f(0:num_k-1))
        f(:)     = 0.0_MK
        
        !------------------------------------------
        ! Calculate all combinations/binomials needed.
        !------------------------------------------

        DO n = 0, num_k-1
           DO s = 1, num_k - 1
              comb_value(n+s,n) = combination(n+s,n)
           END DO
        END DO

        
        !------------------------------------------
        ! coefficients.
        !------------------------------------------
        
        DO k = 0, num_k-1
           
           DO q = 0, k
              
              IF ( p_mark(1,k-1,q) ) THEN
                 f(k) = f(k) + p_value(1,k-q,q) * lambda**q
                 !PRINT *, 2**k*p_value(1,k-q,q)                 
              ELSE
                 f(k) = f(k) + p_func(1,k-q,q) * lambda**q
                 !PRINT *, 2**k*p_func(1,k-q,q)
              END IF
              
              
           END DO
           
           f(k) = f(k) * 2**k
           !PRINT *, f(k)
           !PRINT *, ""
           
        END DO
        
        
        DO p=0, num_p

           s_cur = s_min + p * s_dif
           
           d = 1.0_MK/s_cur/(1.0_MK+lambda)
           
           dd = 1.0_MK
           
           XA11 = 0.0_MK
           XA12 = 0.0_MK
           
           !PRINT *, d, dd, (num_k-1)/2
           DO k = 0, (num_k-1)/2
              
              XA11 = XA11 + f(2*k) * dd
              
              XA12 = XA12 - f(2*k+1) * dd * d
              
              dd = dd*d**2
              
           END DO
           
           
           Fd = 6.0_MK*pi*a*eta*(XA11-XA12)*U
           !Fd = (XA11-XA12)
           
           PRINT *, s_cur-2.0_MK, Fd
        
        END DO
        
9999    CONTINUE
        
        IF ( ASSOCIATED(comb_value)) THEN
           DEALLOCATE(comb_value)
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
        
        IF ( ASSOCIATED(f)) THEN
           DEALLOCATE(f)
        END IF
        
        
      END PROGRAM resistance
       
      
