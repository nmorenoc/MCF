      SUBROUTINE kernel_kernel_Lucy_w(this, rij, w, stat_info)
        !----------------------------------------------------
        !  Subroutine	:  kernel_kernel_Lucy
        !----------------------------------------------------
        !
        !  Purpose      :  Computing Lucy kernel.
        !
        !	 	      	 
        !                
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : V0.1 10.07.2009, original version.
        !
        !----------------------------------------------------
        !  Author       : Xin Bian
        !  Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        
        !----------------------------------------------------
        !  Arguments
        !
        !  rij      : distance of two points.
        !  w        : kernel value.
        ! stat_info : return flag of status.
        !----------------------------------------------------
        
        TYPE(Kernel), INTENT(IN)        :: this
        REAL(MK),  INTENT(IN)           :: rij
        REAL(MK), INTENT(OUT)           :: w    
        INTEGER, INTENT(INOUT)          :: stat_info
        
        
        !----------------------------------------------------
        !  Local variables.
        !----------------------------------------------------

        REAL(MK)                        :: s
        

        !----------------------------------------------------
        !  Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0       
              
        !----------------------------------------------------
        !  Check if the distance is 
        !  non-negative.
        !----------------------------------------------------
        
        IF ( rij < 0 ) THEN
           PRINT *, "kernel_kernel_Lucy : ", &
                "rij should not be negative !"
           stat_info = -1
           GOTO 9999
           
        END IF
        
        
        s = rij / this%h
        
        IF ( s <= 1.0_MK ) THEN
           
           w = this%coef * (1.0_MK+3.0_MK*s)*(1.0_MK-s)**3
           
        ELSE
           
           w = 0.0_MK
           
        END IF
        
9999    CONTINUE        
        
        RETURN        
        
      END SUBROUTINE kernel_kernel_Lucy_w


      SUBROUTINE kernel_kernel_Lucy_w_gradw (this, rij, &
           w, gradW, stat_info)
        !----------------------------------------------------
        ! Subroutine  :  kernel_kerel_Lucy
        !----------------------------------------------------
        !
        ! Purpose     :  Computing w and gradient of
        !                 Lucy kernel.
        !
        !	 	      	 
        !                
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 10.07.2009, original version.
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

        
        !----------------------------------------------------
        !  Arguments
        !
        !  rij      : distance of two points.
        !  w        : kernel value.
        ! gradW     : gradient of kernel.
        ! stat_info : return flag of status.
        !----------------------------------------------------
        
        TYPE(Kernel), INTENT(IN)        :: this
        REAL(MK),  INTENT(IN)           :: rij
        REAL(MK), INTENT(OUT)           :: w    
        REAL(MK), INTENT(OUT)           :: gradW
        INTEGER, INTENT(INOUT)          :: stat_info
        

        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        REAL(MK)                        :: s
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0       
              
        !----------------------------------------------------
        ! Check if the distance is  non-negative.
        !----------------------------------------------------
        
        IF ( rij < 0 ) THEN
           PRINT *, "kernel_kernel_Lucy : ", &
                "rij should not be negative !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        s = rij / this%h
        
        
        IF ( s <= 1.0_MK ) THEN
           
           w     = this%coef * (1.0_MK+3.0_MK*s)*(1.0_MK-s)**3
           gradW = this%coef_grad * 12.0_MK * s * (1.0_MK-s)**2 
           
        ELSE
           
           w     = 0.0_MK
           gradW = 0.0_MK
           
        END IF
        
9999    CONTINUE        
        
        RETURN        
        
      END SUBROUTINE kernel_kernel_Lucy_w_gradw
      
