!----------------------------------------------------
! Functions of ellipse such as r, dr, ddr 
! in polar coordinate system.
! The advantage in polar coordinate:
! it is easy to figure out how an ellipse rotates,
! as the rotated angle appears in the formular
! directly in polar coordinate.
! Reference: wikipedia for basic concepts and
!            for further derivations, check Xin's notes.
!----------------------------------------------------

      REAL(MK) FUNCTION colloid_polar_ellipse_r(a,b,theta,phi)
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        colloid_polar_ellipse_r = &
             SQRT(2.0_MK) * a * b / &
             SQRT( ( b**2-a**2 )* &
             COS(2.0_MK*(theta-phi)) + &
             a**2 + b**2 )
        
      END FUNCTION colloid_polar_ellipse_r
      
      
      REAL(MK) FUNCTION colloid_polar_ellipse_dr(a,b,theta,phi)
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        REAL(MK)                        :: r
        
        r =  colloid_polar_ellipse_r(a,b,theta,phi) 

        colloid_polar_ellipse_dr = &
             r**3.0 * &
             ( b**2-a**2 ) * SIN(2.0_MK*(theta-phi)) / &
             (2.0_MK*a**2*b**2)
        
      END FUNCTION colloid_polar_ellipse_dr
      
      
      REAL(MK) FUNCTION colloid_polar_ellipse_ddr(a,b,theta,phi)
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        REAL(MK)                        :: r, dr
        
        r  =  colloid_polar_ellipse_r(a,b,theta,phi) 
        dr =  colloid_polar_ellipse_dr(a,b,theta,phi)
        
        colloid_polar_ellipse_ddr = &
             3.0_MK * dr** 2 / r  + &
             r**3 * ( b**2-a**2 ) * COS(2.0_MK*(theta-phi)) / &
             a**2/b**2
        
      END FUNCTION colloid_polar_ellipse_ddr
