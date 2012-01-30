#ifdef __GFORTRAN
      SUBROUTINE random_init(this,r_type,d_seed,stat_info)
        !----------------------------------------------------
        ! Subroutine  : random_init
        !----------------------------------------------------
        !
        ! Purpose     : Constructor of Random object,
        !               for initialization of its variables.
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
        INTEGER, INTENT(IN)             :: r_type
        INTEGER, INTENT(IN)             :: d_seed
        INTEGER, INTENT(OUT)            :: stat_info

        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------

        INTEGER, DIMENSION(:),ALLOCATABLE   :: seed
        INTEGER                             :: n, clock,i
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------
        
        stat_info = 0

        this%random_type = r_type
        n                = 1
        
        this%random_uniform_type  = 1
        this%random_Gaussian_type = 1

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
          
        CALL SYSTEM_CLOCK(COUNT=clock)
          
        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)
        
        DEALLOCATE(seed)
        
        RETURN
        
      END SUBROUTINE random_init

#else

      SUBROUTINE random_init(this,r_type,d_seed,stat_info)
        !----------------------------------------------------
        ! Subroutine  : random_init
        !----------------------------------------------------
        !
        ! Purpose     : Constructor of Random object,
        !               for initialization of its variables.
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
        INTEGER, INTENT(IN)             :: r_type
        INTEGER, INTENT(IN)             :: d_seed
        INTEGER, INTENT(OUT)            :: stat_info

        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------

        INTEGER, DIMENSION(:),ALLOCATABLE   :: seed
        INTEGER                         :: k
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------
        
        stat_info = 0

        this%random_type = r_type
        k                = 1
        
        this%random_uniform_type  = 1
        this%random_Gaussian_type = 1

        IF ( d_seed > 0 ) THEN
           CALL RANDOM_SEED(SIZE=k)
           ALLOCATE(seed(k))
           seed(1)  = d_seed
           CALL RANDOM_SEED(PUT=seed(1:k))

        ELSE
           CALL RANDOM_SEED()
           CALL RANDOM_SEED(SIZE=k)
           ALLOCATE(seed(k))
           CALL RANDOM_SEED(GET=seed(1:k))

        END IF

        this%seed = seed(1)
        this%iset = 0

        RETURN
        
      END SUBROUTINE random_init
#endif
      
      SUBROUTINE random_display_parameters(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : random_display_parameters
        !----------------------------------------------------
        !
        ! Purpose     : To display random paramters.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 15.10.2010, original version.
        !
        !----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------

        TYPE(Random),INTENT(IN)         :: this
        INTEGER,INTENT(OUT)             :: stat_info

        
        stat_info = 0
              
        PRINT *, '---***************Start***************---'
        PRINT *, '     Random number parameters'
        PRINT *, '---***********************************---'

        PRINT *, "random type         : ", this%random_type
        PRINT *, "random seed         : ", this%seed
        SELECT CASE ( this%random_type ) 
        CASE (1)
           PRINT *, "random uniform type : ", this%random_uniform_type
        CASE (2)
           PRINT *, "random Gaussian type: ", this%random_Gaussian_type
        END SELECT
        PRINT *, '---****************End****************---'

        
        RETURN
        
      END SUBROUTINE random_display_parameters
