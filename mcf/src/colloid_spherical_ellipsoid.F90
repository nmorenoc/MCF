      REAL(MK) FUNCTION colloid_spherical_ellipsoid_r(a,b,c,theta,phi)
        !------------------------------------------------------------
        ! Purpose   : Function of an ellipsoid, 
        !             distance r on the surface from origin
        !             in spherical coordinate system.
        !------------------------------------------------------------
        
        REAL(MK), INTENT(IN)            :: a
        REAL(MK), INTENT(IN)            :: b
        REAL(MK), INTENT(IN)            :: c
        REAL(MK), INTENT(IN)            :: theta
        REAL(MK), INTENT(IN)            :: phi
        
        colloid_spherical_ellipsoid_r = &
             SQRT( 1.0_MK / ( COS(theta)**2*SIN(phi)**2/a**2 + &
             SIN(theta)**2*SIN(phi)**2/b**2 + &
             COS(phi)**2/c**2 ) )
        
      END FUNCTION colloid_spherical_ellipsoid_r
