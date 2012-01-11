      REAL(MK) FUNCTION random_random(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : random_random
        !----------------------------------------------------
        !
        ! Purpose     : According to the type of random
        !               number generator, return a random
        !               number.
        !
        !
        ! Revisions   : 0.1 03.03. 2009, original version.
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

       
        REAL(MK)                        :: rand
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        !----------------------------------------------------
        ! Select randomber number generator
        !----------------------------------------------------
        
        SELECT CASE(this%random_type)
           
        CASE (1)
           
           DO
              CALL RANDOM_NUMBER(rand)
              
              IF ( ABS(rand) <= 1.0_MK ) THEN
                 EXIT
              END IF
              
           END DO
           
        CASE (2)
           
           DO
              rand = random_random_Gaussian2(this,stat_info_sub)
              
              IF ( ABS(rand) < mcf_random_Gaussian_max_abs ) THEN
                 EXIT
              END IF
              
           END DO
           
           
        END SELECT
        
        random_random = rand
        
        IF(stat_info_sub /=0) THEN
           
           PRINT *, "random_random : ", &
                "random generator has problem !"
           stat_info  = -1
           
        END IF
        
        
        RETURN
        
      END FUNCTION random_random
