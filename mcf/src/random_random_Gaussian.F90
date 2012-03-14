      REAL(MK) FUNCTION random_random_Gaussian1(this,stat_info) 
        !----------------------------------------------------
        ! Subroutine  : random_random_Gaussian1
        !----------------------------------------------------
        !
        ! Purpose     : Return a normal/Gaussian distributed
        !               random number.
        !
        ! Remark      : According to Numerical Recipes 2nd C.
        !
        ! Revisions   : 0.1 06.10. 2010, original version.
        !----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !----------------------------------------------------
        TYPE(Random), INTENT(INOUT)     :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------

        REAL(MK)                        :: v1,v2
        REAL(MK)                        :: fac,rsq
                
        stat_info = 0
        
        IF ( this%iset == 0 ) THEN
           
           DO
           
              CALL RANDOM_NUMBER(v1)
              CALL RANDOM_NUMBER(v2)
              v1 = 2.0_MK * v1 - 1.0_MK
              v2 = 2.0_MK * v2 - 1.0_MK
              rsq = v1**2 + v2**2
              
              !----------------------------------------------
              ! Stop if (v1,v2) are in unit circle.
              !----------------------------------------------

              IF ( rsq > 0.0_MK .AND. rsq < 1.0_MK ) THEN
                 EXIT
              END IF
              
           END DO
           
           fac  = SQRT(-2.0*LOG(rsq)/rsq)
           
           !-------------------------------------------------
           ! Now make the Box-Muller transformation to get
           ! two normal deviates.
           ! Return one and save the other for next time.
           !-------------------------------------------------
           random_random_Gaussian1 = v1 * fac

           this%iset = 1
           this%gset = v2 * fac

        ELSE
           
           this%iset = 0
           random_random_Gaussian1 = this%gset
           
        END IF
        
      END FUNCTION random_random_Gaussian1
      
      
      REAL(MK) FUNCTION random_random_Gaussian2(this,stat_info) 
        !----------------------------------------------------
        ! Subroutine  : random_random_Gaussian1
        !----------------------------------------------------
        !
        ! Purpose     : Return a normal/Gaussian distributed
        !               random number.
        !
        ! Remark      :
        !               The algorithm was found on internet
        !               and can not be found the source again.
        !               However, it proves to work well 
        !               practically.
        !
        ! Revisions   : 0.1 06.10. 2010, original version.
        !----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !----------------------------------------------------
        TYPE(Random), INTENT(IN)        :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------

        REAL(MK), PARAMETER             :: s = 0.449871_MK, &
             t = -0.386595_MK, &
             a = 0.19600_MK, &
             b = 0.25472_MK,&
             r1 = 0.27597_MK, &
             r2 = 0.27846_MK
        REAL(MK)                        ::  u, v, x, y, q
        
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------

        stat_info = this%random_Gaussian_type        
        stat_info = 0
        
        DO
           
           CALL RANDOM_NUMBER(u)
           CALL RANDOM_NUMBER(v)
           
           v = 1.7156_MK * (v - 0.5_MK)
           
           x = u - s
           y = ABS(v) - t
           q = x**2 + y*(a*y - b*x)
           
           !---------------------------------------
           !  Accept P if inside inner ellipse.
           !---------------------------------------
           
           IF (q < r1) THEN
              EXIT
           END IF
           
           !---------------------------------------         
           !  Reject P if outside outer ellipse.
           !---------------------------------------

           IF (q > r2) THEN
              CYCLE
           END IF
           
           !---------------------------------------
           !  Reject P if outside acceptance region.
           !---------------------------------------
           
           IF (v**2 < -4.0_MK*LOG(u)*u**2) THEN
              EXIT
           END IF
           
        END DO
        
        !-------------------------------------
        !  Return ratio of P's coordinates as 
        !  the normal deviate.
        !-------------------------------------
        
        random_random_Gaussian2 = v/u
        
        RETURN
        
      END FUNCTION random_random_Gaussian2
