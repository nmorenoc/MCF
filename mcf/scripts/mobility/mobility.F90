  PROGRAM mobility
    !------------------------------------
    ! Mobility problem calculation
    ! according to Russel, Saville and 
    ! Schowalter 1989<<Colloidal Dispersion>>
    ! or, 
    ! G.K. Batchelor 1976, J. Fluid Mech.(1976), 
    ! vol. 74, part 1 ,pp. 1-29
    !------------------------------------
    
    IMPLICIT NONE
    !------------------------------------
    ! rad : radius of two spheres
    ! h0  : gap between two spheres
    ! r0  : initial distance between centers
    ! rr  : vector connecting two centers
    !------------------------------------
    
    INTEGER, PARAMETER          :: MK= KIND(1.0D0)
    REAL(MK),PARAMETER          :: pi=3.14159265358979_MK
    REAL(MK), DIMENSION(2)      :: rad
    REAL(MK)                    :: h0, r0
    REAL(MK) , DIMENSION(3)     :: rr
    REAL(MK)                    :: lam, rho
    REAL(MK)                    :: density,eta
    REAL(MK) , DIMENSION(2)     :: vol,m
    REAL(MK) , DIMENSION(3,2)   :: fa
    REAL(MK) , DIMENSION(3,2)   :: F
    REAL(MK) , DIMENSION(2,2)   :: A,B
    REAL(MK) , DIMENSION(3,3,2,2):: omega
    REAL(MK)                    :: cof
    INTEGER                     :: i,j,k,p
    REAL(MK) , DIMENSION(3,2)   :: U
    REAL(MK) , DIMENSION(3)     :: U0
    
    rad(1)   = 0.02_MK
    rad(2)   = 0.02_MK
    h0       = 0.01_MK
    r0       = rad(1)+rad(2)+h0
    rr(1)    = r0
    rr(2:3)  = 0.0_MK
    lam      = rad(2)/rad(1)
    rho      = 2*r0/(rad(1)+rad(2))
    vol(1:2) = 4.0*pi*rad(1:2)**3/3.0
    
    density = 1.0e+3_MK
    eta = 0.1_MK
    m(1:2) = vol(1:2) * density
    
    fa(1:3,1:2) = 0.0_MK
    fa(1,1) = 9.375e-4_MK;
    fa(1,2) = fa(1,1)
    
    F(1:3,1) = m(1)* fa(1:3,1)
    F(1:3,2) = m(2)* fa(1:3,2)
    U0(1:3)  = F(1:3,1) / (6.0*pi*eta*rad(1))
    
    A(1,1) = 1.0_MK- 60.0_MK*lam**3/(1+lam)**4/rho**4
    A(1,2) = 3.0_MK/2.0_MK/rho - &
         2.0_MK*(1.0_MK+lam**2)/(1.0_MK+lam)**2/rho**3
    A(2,1) = A(1,2)
    A(2,2) = A(1,1)
    
    B(1,1) = 1
    B(1,2) = 3.0_MK/4.0_MK/rho + (1.0_MK+lam**2) / (1+lam)**2/rho**3
    B(2,1) = B(1,2)
    B(2,2) = B(1,1)
    
    omega(:,:,:,:) = 0.0_MK
    cof = 3.0*pi*eta*(rad(1)+rad(2))
    
    DO j = 1, 2
       DO i = 1, 2
          
          DO p = 1, 3
             DO k =1, 3
                
                omega(k,p,i,j) = &
                     ( A(i,j)-B(i,j)) * rr(k) * rr(p) / r0**2
                
                IF (k==p)  THEN
                   omega(k,p,i,j) = omega(k,p,i,j) + B(i,j)
                ENd IF
                
                omega(k,p,i,j)= omega(k,p,i,j) / cof
                
             END DO
          END DO
          
       END DO
    END DO
    
    U(1:3,1:2) = 0.0_MK
    
    !------------------------------------
    ! Print out data onto screen
    !------------------------------------
    
    PRINT *, "rr"
    PRINT *, rr(1:3)
    
    PRINT *, "Coefficients of mobility tensor"
    PRINT *, "A : "
    PRINT *, A(1,1), A(1,2)
    PRINT *, A(2,1), A(2,2)

    PRINT *, "B : "
    PRINT *, B(1,1), B(1,2)
    PRINT *, B(2,1), B(2,2)
    
    
    DO i= 1,2
       DO j =1,2
          U(1:3,i) = U(1:3,i) + &
               MATMUL(omega(1:3,1:3,i,j),F(1:3,j))
       END DO
    END DO
    
    PRINT *, "Mobility tensor"
    DO j = 1,2
       DO i =1,2
          PRINT *, "omega(", i,j,") :"
          PRINT *, omega(1:3,1:3,i,j)
       END DO
    END DO
    
    PRINT *, "Force on two spheres"
    
    PRINT *, "F : "
    PRINT *, F(1,1), F(1,2)
    PRINT *, F(2,1), F(2,2)
    PRINT *, F(3,1), F(3,2)
    
    PRINT *, "Velocity of two spheres"
    
    PRINT *, "U : "
    PRINT *, U(1,1), U(1,2)
    PRINT *, U(2,1), U(2,2)
    PRINT *, U(3,1), U(3,2)
  
    PRINT *, "Settling velocity of isolate sphere"
    PRINT *, "U0: "
    PRINT *, U0(1)
    PRINT *, U0(2)
    PRINT *, U0(3)

    
  END PROGRAM mobility
  
  
