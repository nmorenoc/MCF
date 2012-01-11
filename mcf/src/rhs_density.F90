      SUBROUTINE  rhs_density_ff(this,w,mi,mj,rhoi,rhoj,stat_info)
        !----------------------------------------------------
        ! Fluid Fluid interaction for density summation.
        !----------------------------------------------------

        TYPE(Rhs), INTENT(IN)           :: this
        REAL(MK), INTENT(IN)            :: w
        REAL(MK), INTENT(IN)            :: mi
        REAL(MK), INTENT(IN)            :: mj
        REAL(MK), INTENT(OUT)           :: rhoi
        REAL(MK), INTENT(OUT)           :: rhoj
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        INTEGER                         :: stat_info_sub
        
        
        stat_info     = 0
        stat_info_sub = 0

        rhoi = 0.0_MK
        rhoj = 0.0_MK
        
        SELECT CASE(this%rhs_density_type)

        CASE (1)
           
           CALL rhs_density_ff_phys(this,w,mi,mj,rhoi,rhoj,stat_info)
           
        CASE (2)
           
           CALL rhs_density_ff_num(this,w,rhoi,rhoj,stat_info)
           
        END SELECT

        IF (stat_info_sub /= 0) THEN
           PRINT *, "rhs_density_ff : ", &
                "Density calculation has problem !"
           stat_info = -1
           GOTO 9999

        END IF

9999    CONTINUE
        
        
        RETURN 

      END SUBROUTINE rhs_density_ff
      
      
      SUBROUTINE  rhs_density_ff_phys(this,w,mi,mj,rhoi,rhoj,stat_info)

        TYPE(Rhs), INTENT(IN)           :: this
        REAL(MK), INTENT(IN)            :: w
        REAL(MK), INTENT(IN)            :: mi
        REAL(MK), INTENT(IN)            :: mj
        REAL(MK), INTENT(INOUT)         :: rhoi
        REAL(MK), INTENT(INOUT)         :: rhoj
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF ( this%rhs_density_type /= 1 ) THEN
           stat_info = -1
           GOTO 9999
        END IF

        rhoi = mj * w
        rhoj = mi * w
        
9999    CONTINUE
        
        RETURN 
        
      END SUBROUTINE rhs_density_ff_phys
      

      SUBROUTINE  rhs_density_ff_num(this,w,numi,numj,stat_info)
        
        TYPE(Rhs), INTENT(IN)           :: this
        REAL(MK), INTENT(IN)            :: w
        REAL(MK), INTENT(INOUT)         :: numi
        REAL(MK), INTENT(INOUT)         :: numj
        
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0

        IF ( this%rhs_density_type /= 2 ) THEN
           stat_info = -1
           GOTO 9999
        END IF

        numi = w
        numj = w
        
9999    CONTINUE
        
        RETURN 
        
      END SUBROUTINE rhs_density_ff_num
      
