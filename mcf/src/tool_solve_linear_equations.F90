      SUBROUTINE tool_solve_linear_equations(this, &
           N, NRHS, A, LDA, IPIV, B, LDB,stat_info)
        !----------------------------------------------------
        ! Subroutine  : tool_solve_linear_equations
        !----------------------------------------------------
        !
        ! Purpose     : Solve a system of linear equations
        !
        ! Routines    :
        !
        ! Remarks     : Check Lapack for remark of parameters
        !
        ! References  : Lapack 
        !
        ! Revisions   : V0.1 30.04 2012, original version.
        !
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments     
        !----------------------------------------------------
        
        TYPE(Tool), INTENT(IN)                  :: this
        INTEGER, INTENT(IN)                     :: N
        INTEGER, INTENT(IN)                     :: NRHS
        REAL(MK),DIMENSION(:,:), INTENT(IN)     :: A
        INTEGER, INTENT(IN)                     :: LDA
        INTEGER, DIMENSION(:), INTENT(INOUT)    :: IPIV
        REAL(MK),DIMENSION(:,:), INTENT(INOUT)  :: B
        INTEGER, INTENT(IN)                     :: LDB
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization
        !
        ! This is supposed to be used, otherwise,
        ! compiler complains that it is not used.
        !----------------------------------------------------
        
        stat_info     = this%flag
        stat_info     = 0
        stat_info_sub = 0
        
        stat_info = -1
        PRINT *, __FILE__,__LINE__,&
             "Lapack is not linked!"
        GOTO 9999
        
        IF ( MK == ppm_kind_single ) THEN
           
           !CALL SGESV(N,NRHS,A,LDA,IPIV,B,LDB,stat_info_sub)
           
        ELSE IF ( MK == ppm_kind_double ) THEN
           
           !CALL DGESV(N,NRHS,A,LDA,IPIV,B,LDB,stat_info_sub)
           
        END IF

        IF ( stat_info_sub /= 0 ) THEN
           
           PRINT *, "N, LDA, LDB: ", N, LDA, LDB
           PRINT *, "A, B:", A(:,:), B(:,:)
           PRINT *, __FILE__, __LINE__, "does not work!"
           stat_info = -1
           GOTO 9999
           
        END IF
        
        
        !----------------------------------------------------
        ! Return.
        !----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE tool_solve_linear_equations
      
      
