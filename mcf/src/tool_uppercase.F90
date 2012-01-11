      SUBROUTINE tool_uppercase(this,string,ilen,stat_info)
        !----------------------------------------------------
        !  Subroutine   :     tool_uppercase
        !----------------------------------------------------
        !
        !  Purpose      : Converts a character string to 
        !                 upper case letters.
        !
        !  Routines     :
        !
        !  Remarks      :
        !
        !  References   :
        !
        !  Revisions    : V0.1 03.03 2009 original version.
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
        
        TYPE(Tool), INTENT(IN)          :: this
        CHARACTER(LEN=*), INTENT(INOUT) :: string
        INTEGER, INTENT(IN)             :: ilen
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables 
        !----------------------------------------------------
        
        INTEGER               :: i,j
        INTEGER               :: i1,i2,i3,iadd
        
        
        !----------------------------------------------------
        ! Initialise
        ! this is supposed to be used, otherwise,
        ! compiler complains that it is not used.
        !----------------------------------------------------
        
        stat_info = this%flag
        stat_info = 0
        
        !----------------------------------------------------
        ! Convert to upper case
        !----------------------------------------------------
        
        i1   = IACHAR('a') - 1
        i2   = IACHAR('z') + 1
        i3   = IACHAR('A')
        iadd = i3 - i1 - 1
        
        DO i=1,ilen
           j = IACHAR(string(i:i))
           IF (j.GT.i1.AND.j.LT.i2) THEN
              string(i:i) = CHAR(j+iadd)
           END IF
        END DO
        
        !----------------------------------------------------
        ! Return 
        !----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE tool_uppercase
      
