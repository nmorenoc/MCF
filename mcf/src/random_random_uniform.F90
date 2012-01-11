      REAL(MK) FUNCTION random_random_uniform(this,stat_info) 
        !----------------------------------------------------
        ! Subroutine  : random_random_uniform
        !----------------------------------------------------
        !
        ! Purpose     : Return a uniform distributed
        !               random number.
        !
        ! Remark      : 
        !
        ! Revisions   : 0.1 19.10. 2010, original version.
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
        
        INTEGER                         :: stat_info_sub
        REAL(MK)                        :: rand
        INTEGER, PARAMETER :: K4B=selected_int_kind(9)
        INTEGER(K4B)                    :: idum

        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        rand =  -1.0_MK
        
        SELECT CASE(this%random_uniform_type)

        CASE (1)
           
           rand = random_random_uniform1(this,stat_info_sub)
           
        CASE (2)
           
           idum = -1
           rand = random_random_uniform2(this,idum,stat_info_sub)

        CASE DEFAULT
           
           PRINT *, "random_random_uniform : type not know "
           stat_info = -1
           GOTO 9999
           
        END SELECT
        
9999    CONTINUE
        
        random_random_uniform = rand
        
      END FUNCTION random_random_uniform


      REAL(MK) FUNCTION random_random_uniform1(this,stat_info) 
        !----------------------------------------------------
        ! Subroutine  : random_random_uniform1
        !----------------------------------------------------
        !
        ! Purpose     : Return a uniform distributed
        !               random number.
        !
        ! Remark      : Use Fortran 90 intrinct function..
        !
        ! Revisions   : 0.1 19.10. 2010, original version.
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
        
        REAL(MK)                        :: rand
        
        stat_info = 0
        
        CALL RANDOM_NUMBER(rand)
        
        random_random_uniform1 = rand
        
      END FUNCTION random_random_uniform1
      
      
      REAL(MK) FUNCTION random_random_uniform2(this,idum,stat_info) 
        !----------------------------------------------------
        ! Subroutine  : random_random_uniform2
        !----------------------------------------------------
        !
        ! Purpose     : Return a uniform distributed
        !               random number.
        !
        ! Remark      : According to Numerical Recipes 2nd 
        !               Fortran 90.
        !
        ! Revisions   : 0.1 19.10. 2010, original version.
        !----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        IMPLICIT NONE
        INTEGER, PARAMETER :: K4B=selected_int_kind(9)

        !----------------------------------------------------
        ! Arguments
        !----------------------------------------------------
        
        TYPE(Random), INTENT(IN)        :: this
	INTEGER(K4B), INTENT(INOUT)     :: idum
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER(K4B), PARAMETER         :: IA=16807,IM=2147483647
        INTEGER(K4B), PARAMETER         :: IQ=127773,IR=2836

        REAL(MK), SAVE                  :: am
        INTEGER(K4B), SAVE              :: ix=-1,iy=-1,k
        
        !----------------------------------------------------
        ! Initialization.
        !----------------------------------------------------

        stat_info = 0
        
        if (idum <= 0 .or. iy < 0) then
           am=nearest(1.0,-1.0)/IM
           iy=ior(ieor(888889999,abs(idum)),1)
           ix=ieor(777755555,abs(idum))
           idum=abs(idum)+1
        end if
        ix=ieor(ix,ishft(ix,13))
        ix=ieor(ix,ishft(ix,-17))
        ix=ieor(ix,ishft(ix,5))
        k=iy/IQ
        iy=IA*(iy-k*IQ)-IR*k
        if (iy < 0) iy=iy+IM
        random_random_uniform2=am*ior(iand(IM,ieor(ix,iy)),1)
       
      END FUNCTION random_random_uniform2
