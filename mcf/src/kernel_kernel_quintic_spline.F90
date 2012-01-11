      SUBROUTINE kernel_kernel_quintic_spline_w(this, rij, &
           w, stat_info)
        !----------------------------------------------------
        ! Subroutine  :  kernel_kerel_quintic_spline_w
        !----------------------------------------------------
        !
        ! Purpose     :  Computing quintic spline kernel.
        !
        !	 	      	 
        !                
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 10.07.2009, rewritten version.
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
        ! Arguments
        !
        ! rij       : distance of two points.
        ! w         : kernel value.
        ! stat_info : return flag of status.
        !----------------------------------------------------
        
        TYPE(Kernel), INTENT(IN)        :: this
        REAL(MK),  INTENT(IN)           :: rij
        REAL(MK), INTENT(OUT)           :: w    
        INTEGER, INTENT(INOUT)          :: stat_info
        
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        REAL(MK)                        :: s
        REAL(MK)                        :: s1,s2,s3
        REAL(MK)                        :: s1_4,s2_4,s3_4
        REAL(MK)                        :: s1_5,s2_5,s3_5

        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0       
        
        !----------------------------------------------------
        ! Check if the distance is non-negative.
        !----------------------------------------------------
        
        IF ( rij < 0 ) THEN
           PRINT *, "kernel_kernel_quintic_spline : ", &
                "rij should not be negative !"
           stat_info = -1
           GOTO 9999
        END IF
        
        s    = rij / this%h
        s1   = 1.0_MK - s
        s2   = 2.0_MK - s
        s3   = 3.0_MK - s
        
        s1_4  = s1**4
        s2_4  = s2**4
        s3_4  = s3**4
        
        s1_5  = 15.0_MK * s1_4 * s1
        s2_5  = -6.0_MK * s2_4 * s2
        s3_5  = s3_4 * s3
        
        
        IF ( s < 1.0_MK) THEN
           
           w = this%coef * (s3_5  + s2_5 +  s1_5)
           
        ELSE IF ( s < 2.0_MK ) THEN

           w = this%coef * (s3_5 + s2_5 )
           
        ELSE IF( s < 3.0_MK ) THEN
           
           w = this%coef * s3_5
           
        ELSE
           
           w = 0.0_MK
           
        END IF
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE kernel_kernel_quintic_spline_w
      
      
      SUBROUTINE kernel_kernel_quintic_spline_w_gradw(this, rij, &
           w, gradW, stat_info)
        !----------------------------------------------------
        ! Subroutine  :  kernel_kerel_quintic_spline
        !----------------------------------------------------
        !
        ! Purpose     :  Computing dirivative of 
        !                quintic spline kernel.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 10.07.2009, rewritten version.
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
        ! Arguments
        !
        ! rij      : distance of two points.
        ! w        : kernel value.
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
        REAL(MK)                        :: s1,s2,s3
        REAL(MK)                        :: s1_4,s2_4,s3_4
        REAL(MK)                        :: s1_5,s2_5,s3_5

        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0       
        
        !----------------------------------------------------
        ! Check if the distance is non-negative.
        !----------------------------------------------------
        
        IF ( rij < 0 ) THEN
           PRINT *, "kernel_kernel_quintic_spline : ", &
                "rij should not be negative !"
           stat_info = -1
           GOTO 9999
        END IF
        
        s    = rij / this%h
        s1   = 1.0_MK - s
        s2   = 2.0_MK - s
        s3   = 3.0_MK - s
        
        s1_4  = s1**4
        s2_4  = s2**4
        s3_4  = s3**4
        
        s1_5  = 15.0_MK * s1_4 * s1
        s2_5  = -6.0_MK * s2_4 * s2
        s3_5  = s3_4 * s3
        
        s1_4  =  75.0_MK * s1_4
        s2_4  = -30.0_MK * s2_4
        s3_4  =  5.0_MK  * s3_4
        
        
        IF ( s < 1.0_MK) THEN
           
           w     = this%coef * (s3_5  + s2_5 +  s1_5)
           gradW = this%coef_grad * (s3_4 + s2_4 + s1_4)
           
        ELSE IF ( s < 2.0_MK ) THEN

           w     = this%coef * (s3_5 + s2_5 )
           gradW = this%coef_grad * (s3_4 + s2_4 )
           
        ELSE IF( s < 3.0_MK ) THEN
           
           w     = this%coef * s3_5
           gradW = this%coef_grad * s3_4
           
        ELSE
           
           w     = 0.0_MK
           gradW = 0.0_MK
           
        END IF
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE kernel_kernel_quintic_spline_w_gradw
