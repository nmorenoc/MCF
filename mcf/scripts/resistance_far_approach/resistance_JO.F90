      PROGRAM resisance
        !----------------------------------------------------
        ! Resistance problem calculation
        ! according to Jeffrey and Onishi 1984, or
        ! Kim and Karrila <<microdydrodynamics>> 1991.
        !----------------------------------------------------
        
        IMPLICIT NONE
        
        !----------------------------------------------------
        ! a1,a2,a : radius of two spheres
        ! R       : initial distance between centers
        ! b       : a2/a1
        !----------------------------------------------------
        
        INTEGER, PARAMETER          :: MK= KIND(1.0D0)
        REAL(MK),PARAMETER              :: pi=3.14159265358979_MK
        REAL(MK)                        :: a
        REAL(MK),DIMENSION(0:10)        :: b
        REAL(MK)                        :: R
        REAL(MK)                        :: d,dd
        REAL(MK),DIMENSION(0:11)        :: f
        INTEGER                         :: k
        REAL(MK)                        :: eta
        REAL(MK)                        :: XA11
        REAL(MK)                        :: XA12
        REAL(MK)                        :: U
        REAL(MK)                        :: Fd
        CHARACTER(256)                  :: cbuf
        INTEGER                         :: n_arg
        

        
        R = 0.08_MK
        a = 0.02_MK
        
        eta  = 0.1_MK
        XA11 = 0.0_MK
        XA12 = 0.0_MK
        
        U    = 1.25e-4_MK
        Fd   = 0.0_MK
        
        n_arg = IARGC()
        
        IF ( n_arg > 0 ) THEN
           CALL GETARG(1,cbuf)
           READ(cbuf,*), R
        END IF
        
        d = a / ( 2.0_MK * R )

        !-----------------------
        ! power of b.
        !-----------------------
        
        b(0) = 1.0_MK
        b(1) = 1.0_MK
        
        DO k = 2, 10
           b(k) = b(k-1)*b(1)
        END DO
        
        
        !-----------------------
        ! coefficients.
        !-----------------------

        f(0) = b(0)
        f(1) = 3.0_MK*b(1)
        f(2) = 9.0_MK*b(1)
        f(3) = -4.0_MK*b(1) + 27.0_MK*b(2) -4.0_MK*b(3)
        f(4) = -24.0_MK*b(1) + 81.0_MK*b(2)+36.0_MK*b(3)
        f(5) = 72.0_MK*b(2) + 243.0_MK*b(3) + 72_MK*b(4)
        f(6) = 16.0_MK*b(1) + 108.0_MK*b(2)+281.0_MK*b(3) + &
             648.0_MK*b(4) + 144.0_MK*b(5)
        f(7) = 288.0_MK*b(2)+ 1620.0_MK*b(3) + &
             1515.0_MK*b(4) + 1620.0_MK*b(5) + 288.0_MK*b(6)
        f(8) = 576.0_MK*b(2)+ 4848.0_MK*b(3) + &
             5409.0_MK*b(4) + 4524.0_MK*b(5) + 3888.0_MK*b(6) + &
             576.0_MK*b(7)
        f(9) = 1152.0_MK*b(2)+ 9072.0_MK*b(3) + &
             14752.0_MK*b(4) + 26163.0_MK*b(5) + 14752.0_MK*b(6) + &
             9072.0_MK*b(7) + 1152.0_MK*b(8)
        f(10) = 2304.0_MK*b(2)+ 20736.0_MK*b(3) + &
             42804.0_MK*b(4) + 115849.0_MK*b(5) + 76176.0_MK*b(6) + &
             39264.0_MK*b(7) + 20736.0_MK*b(8) + 2304.0_MK*b(9)
        f(11) = 4608.0_MK*b(2)+ 46656.0_MK*b(3) + &
             108912.0_MK*b(4) + 269100.0_MK*b(5) + 319899.0_MK*b(6) + &
             269100.0_MK*b(7) + 108912.0_MK*b(8) + 46656.0_MK*b(9) + &
             4608.0_MK*b(10)
        
        dd = 1.0_MK
        
        DO k = 0, 5
           
           XA11 = XA11 + f(2*k) * dd
           
           XA12 = XA12 + f(2*k+1) * dd * d
           
           dd = dd*d**2
           
        END DO

        XA11 = 6.0_MK * pi * a * XA11
        XA12 = -6.0_MK * pi * a * XA12

        Fd = eta * (XA11-XA12) * U

        PRINT *, "XA11   : ", XA11
        PRINT *, "XA12   : ", XA12
        PRINT *, "F drag : ", Fd
        
        
9999    CONTINUE
        
        
      END PROGRAM resisance
      
      
