      PROGRAM resisance
        !----------------------------------------------------
        ! Resistance problem calculation
        ! according to Russel, Saville and 
        ! Schowalter 1989<<Colloidal Dispersion>>, pp 45-51.
        ! or 
        ! G.K. Batchelor 1976, J. Fluid Mech.(1976), 
        ! vol. 74, part 1 ,pp. 1-29
        !----------------------------------------------------
        
        IMPLICIT NONE
        
        !----------------------------------------------------
        ! rad : radius of two spheres
        ! h0  : gap between two spheres
        ! r0  : initial distance between centers
        ! rr  : vector connecting two centers
        !----------------------------------------------------
        
        INTEGER, PARAMETER          :: MK= KIND(1.0D0)
        REAL(MK),PARAMETER          :: pi=3.14159265358979_MK
        REAL(MK), DIMENSION(2)      :: rad
        REAL(MK)                    :: h0, r0
        REAL(MK) , DIMENSION(3)     :: rr
        REAL(MK)                    :: lam, rho
        REAL(MK)                    :: eta
        REAL(MK) , DIMENSION(3,2)   :: U
        REAL(MK) , DIMENSION(3,2)   :: F
        REAL(MK) , DIMENSION(2,2)   :: A,B
        REAL(MK) , DIMENSION(3,3,2,2):: omega
        REAL(MK) , DIMENSION(3,3,2,2):: fc
        REAL(MK) , DIMENSION(3,3)    :: mat1,mat2
        
        REAL(MK)                    :: cof
        INTEGER                     :: i,j,k,p
        REAL(MK) , DIMENSION(3)     :: F0
        
        CHARACTER(256)                  :: cbuf
        INTEGER                         :: n_arg
        
        INTEGER                     :: stat_info
        
        stat_info = 0
        
        rad(1)   = 0.02_MK
        rad(2)   = 0.02_MK
        h0       = 0.04_MK
        
        n_arg = IARGC()

        IF ( n_arg > 0 ) THEN
           CALL GETARG(1,cbuf)
           READ(cbuf,*), h0
           PRINT *, "h0 : ", h0
        END IF
        
        r0       = rad(1)+rad(2)+h0
        rr(1)    = r0
        rr(2:3)  = 0.0_MK
        lam      = rad(2)/rad(1)
        rho      = 2*r0/(rad(1)+rad(2))
        
        
       
        !----------------------------------------------------
        ! Print data onto screen.
        !----------------------------------------------------
        
        PRINT *, "lamda : ", lam
        PRINT *, "rho   : ", rho

        U(1:3,1:2) = 0.0_MK
        U(1,1)     = 1.25e-4_MK
        U(1,2)     = -1.25e-4_MK
        
        eta      = 0.1
        F0(1:3)  = -6.0 * pi * eta*rad(1) * U(1:3,1)
        
        A(1,1) = 1.0_MK- 60.0_MK*lam**3/(1+lam)**4/rho**4
        A(1,2) = 3.0_MK/2.0_MK/rho - &
             2.0_MK*(1.0_MK+lam**2)/(1.0_MK+lam)**2/rho**3
        A(2,1) = A(1,2)
        A(2,2) = A(1,1)
        
        B(1,1) = 1
        B(1,2) = 3.0_MK/4.0_MK/rho + (1.0_MK+lam**2) / (1+lam)**2/rho**3
        B(2,1) = B(1,2)
        B(2,2) = B(1,1)
        
        !----------------------------------------------------
        ! Print data onto screen.
        !----------------------------------------------------
        
        PRINT *, "Coefficients of mobility tensor"
        PRINT *, "A : "
        PRINT *, A(1,1), A(1,2)
        PRINT *, A(2,1), A(2,2)
        
        PRINT *, "B : "
        PRINT *, B(1,1), B(1,2)
        PRINT *, B(2,1), B(2,2)
        
        
        !----------------------------------------------------
        ! Caculating mobility tensor.
        !----------------------------------------------------
        
        omega(:,:,:,:) = 0.0_MK
        cof = 3.0_MK*pi*eta*(rad(1)+rad(2))
        
        DO j = 1, 2
           DO i = 1, 2
              
              DO p = 1, 3
                 DO k =1, 3
                    
                    omega(k,p,i,j) = &
                     ( A(i,j)-B(i,j)) * rr(k) * rr(p) / r0**2
                    
                    IF ( k==p )  THEN
                       
                       omega(k,p,i,j) = omega(k,p,i,j) + B(i,j)
                       
                    END IF
                    
                    omega(k,p,i,j)= omega(k,p,i,j) / cof
                    
                 END DO
              END DO
              
           END DO
           
        END DO
        
        !----------------------------------------------------
        ! Print data onto screen.
        !----------------------------------------------------
        
        PRINT *, "Mobility tensor"
        
        DO j = 1,2
           DO i =1,2
              PRINT *, "omega(", i,j,") :"
              PRINT *, omega(1:3,1:3,i,j)
           END DO
        END DO
        
        !----------------------------------------------------
        ! Calculating resistance tensor.
        !----------------------------------------------------
        
        fc(:,:,:,:) = 0.0_MK
        
        DO i = 1, 2
           
           DO j = 1, 2
              
              IF( i == j ) THEN
                 
                 p = i
                 k = 3-i
                 
                 CALL  matrix_inverse(omega(1:3,1:3,k,k), &
                      mat1(1:3,1:3),3,stat_info)
                 
                 IF ( stat_info /=0 ) THEN
                    PRINT *, "inversing matrix has problem !"
                    GOTO 9999
                 END IF
                 
                 mat2(1:3,1:3) = omega(1:3,1:3,p,p) - &
                      MATMUL(MATMUL(omega(1:3,1:3,p,k), &
                      mat1(1:3,1:3)), omega(1:3,1:3,k,p))
                 
              ELSE
                 
                 p = i
                 k = j
                 
                 CALL  matrix_inverse(omega(1:3,1:3,p,k), &
                      mat1(1:3,1:3),3,stat_info)
                 
                 mat2(1:3,1:3) = omega(1:3,1:3,k,p) - &
                      MATMUL(MATMUL(omega(1:3,1:3,k,k), &
                      mat1(1:3,1:3)), omega(1:3,1:3,p,p))
                 
              END IF ! i == j
              
              CALL matrix_inverse(mat2(1:3,1:3), &
                   fc(1:3,1:3,i,j),3,stat_info)
              
              IF ( stat_info /=0 ) THEN
                 PRINT *, "inversing matrix has problem !"
                 GOTO 9999
              END IF
              
           END DO ! j
           
        END DO ! i
        
        !----------------------------------------------------
        ! Print data onto screen.
        !----------------------------------------------------
        
        
        PRINT *, "Friction/Resistance tensor"
        
        DO j = 1,2
           DO i =1,2
              PRINT *, "fc(", i,j,") :"
              PRINT *, fc(1:3,1:3,i,j)
           END DO
        END DO
        
        
        !----------------------------------------------------
        ! Print data onto screen.
        !----------------------------------------------------
        
        PRINT *, "Velocity of two spheres"
        
        PRINT *, "U : "
        PRINT *, U(1,1), U(1,2)
        PRINT *, U(2,1), U(2,2)
        PRINT *, U(3,1), U(3,2)
        
         
        PRINT *, "Settling drag of an isolated sphere"
        
        PRINT *, "F0: "
        PRINT *, F0(1)
        PRINT *, F0(2)
        PRINT *, F0(3)

        !----------------------------------------------------
        ! Calculate resistance/drag.
        !----------------------------------------------------
        
        F(1:3,1:2) = 0.0_MK
        
        DO i= 1,2
           DO j =1,2
              F(1:3,i) = F(1:3,i) + &
                   MATMUL(fc(1:3,1:3,i,j),U(1:3,j))
           END DO
        END DO
        
        
        PRINT *, "Force on two spheres"
        
        PRINT *, "F : "
        PRINT *, -F(1,1), -F(1,2)
        PRINT *, -F(2,1), -F(2,2)
        PRINT *, -F(3,1), -F(3,2)
        
9999    CONTINUE
        
      END PROGRAM resisance
      
      
      
      SUBROUTINE matrix_inverse(matrix, inverse, n, errorflag)
        
        IMPLICIT NONE
        
        !------------------------------------
        ! Inverse diagonal matrix
        !------------------------------------
        
        ! Declarations
        INTEGER, PARAMETER                   :: MK= KIND(1.0D0)
        
        REAL(MK),INTENT(IN), DIMENSION(n,n)  :: matrix  !Input matrix
        REAL(MK),INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix
        INTEGER, INTENT(IN)                  :: n
        
        !Return error status. -1 for error, 0 for normal
        INTEGER, INTENT(OUT)                :: errorflag  
        
        INTEGER                             :: i,j
        
        errorflag   = 0
        inverse(:,:)= 0.0_MK
        
        DO j=1,n
           
           DO i=1,n
              
              IF ( i==j ) THEN
                 
                 IF ( matrix(i,j) == 0.0_MK ) THEN
                    
                    PRINT *, "Diagonal element is zero !"
                    errorflag = -1
                    GOTO 9999
                    
                 ELSE
                    
                    inverse(i,j) = 1.0_MK/matrix(i,j)
                    
                 END IF
                 
              END IF
              
           END DO
           
        END DO
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE matrix_inverse
      
      
