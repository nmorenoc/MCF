      SUBROUTINE colloid_compute_repulsion_cw(this,&
           x,sid,F,FB,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_repulsion_cw
        !----------------------------------------------------
        !
        ! Purpose     : The gap between near contacting
        !               colloid and wall can be too small
        !               and they tend to overlap.
        !               Compute an extra repulsive force
        !               which prevents them to overlap
        !               or become too close.
        !
        ! Routines    :
        !
        ! References  : 1) Brady and Bossis,
        !                  J. Fluid Mech. vol. 155,
        !                  pp. 105-129.1985
        !               2) Ball and Melrose, 
        !                  Adv. Colloid Interface Sci.
        !                  59 19-30, 1995.
        !               3) Dratler and Schowalter,
        !                  J. Fluid Mech.
        !                  vol. 325, pp 53-77. 1996.
        !               4) Sierou and Brady,
        !                  J. Rheol. 46(5), 1031-1056, 2002.
        !              
        !
        ! Remarks     : V0.3 6.3.2012, change the second
        !               repulsive force, assuming radius
        !               of colloid is universally one.
        !
        !               V0.2 19.11.2010, one pair version.
        !
        !               V0.1 5.11 2010, original version,
        !               for all pairs.
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
        ! Arguments:
        !
        ! input:
        !   x: location of colloid
        ! sid: species id of colloid
        !
        ! output:
        ! F : force on colloid
        ! FB: force on wall boundary
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:),INTENT(IN)       :: x
        INTEGER, INTENT(IN)                     :: sid
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: F
        REAL(MK), DIMENSION(:,:),INTENT(OUT)    :: FB
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !
        ! s_drag   : drag on this process
        ! t_drag   : total drag on all processes
        ! s_torque : torque on this process
        ! t_torque : total torque on all processes        
        !----------------------------------------------------
        
        INTEGER                                  :: dim,num
        INTEGER                                  :: i,j
        REAL(MK)                                 :: hn,hm,h1,h2
        REAL(MK)                                 :: F0,Ft
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        
        dim = this%num_dim
        num = dim * 2
        
        IF ( SIZE(F,1) /= dim ) THEN
           PRINT *, "colloid_compute_repulsion_cw : ", &
                "input F dimension does not match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( SIZE(FB,1) /= dim ) THEN
           PRINT *, "colloid_compute_repulsion_cw : ", &
                "input FB dimension does not match !"
           stat_info = -1
           GOTO 9999
        END IF

        IF ( SIZE(FB,2) /= num ) THEN
           PRINT *, "colloid_compute_repulsion_cw : ", &
                "input FB number does not match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
        F(:)    = 0.0_MK
        FB(:,:) = 0.0_MK
        
        hn = this%cw_repul_cut_off
        hm = this%cw_repul_cut_on
        F0 = this%cw_repul_F0
        Ft = 0.0_MK
        
        !----------------------------------------------------
        ! Loop each dimension of solid wall.
        !----------------------------------------------------
        
        DO j = 1, dim
           
           IF ( this%bcdef(2*j-1) == ppm_param_bcdef_wall_solid ) THEN
              
              !---------------------------------------------
              ! Calculate distance of the colloid to
              ! each wall.
              !----------------------------------------------
              
              h1 = x(j)-this%min_phys(j)-this%radius(1,sid)
              h2 = this%max_phys(j)-x(j)-this%radius(1,sid)
              
              !----------------------------------------------
              ! distance to down wall smaller than 5*hn
              !----------------------------------------------
              
              IF ( h1 < 5.0_MK * hn ) THEN
                 
                 SELECT CASE (this%cw_repul_type)
                    
                 CASE (mcf_cw_repul_type_Hookean)
                    !----------------------------------------
                    ! For linear spring force, it is
                    ! zero at hn and beyond.
                    !----------------------------------------
                    
                    IF ( h1 < hn ) THEN
                       
                       !-------------------------------------
                       ! If gap smaller than minimal 
                       ! allowed gap, set it to minimum.
                       !-------------------------------------
                       
                       IF ( h1 < hm ) THEN
                          
                          h1 = hm
                          
                       END IF
                       
                       Ft = F0 - F0*h1/hn
                       
                    END IF
                    
                 CASE (mcf_cw_repul_type_DLVO)
                    
                    !----------------------------------------
                    ! For DLVO force, it is not zero at hn
                    ! and less than 1/100 beyond 5*hn.
                    !----------------------------------------
                       
                    !----------------------------------------
                    ! If gap smaller than minimal 
                    ! allowed gap, set it to minimum.
                    !----------------------------------------
                    
                    IF ( h1 < hm ) THEN
                       
                       h1 = hm
                       
                    END IF
                    
                    Ft = F0 /hn * EXP(-h1/hn) /(1.0_MK-EXP(-h1/hn))
                    
                 END SELECT ! cw_repul_type
                 
                 F(j) = F(j) +  Ft
                 
                 !-------------------------------------------
                 ! Collect force on boundary
                 !-------------------------------------------
                 
                 FB(j,2*j-1) = -Ft
                 
              END IF ! h1 < 5.0*hn
              
              
              IF ( h2 < 5.0_MK * hn ) THEN
                 
                 SELECT CASE (this%cw_repul_type)
                    
                 CASE (mcf_cw_repul_type_Hookean)
                    
                    !----------------------------------------
                    ! For linear spring force, it is
                    ! zero at hn and beyond.
                    !----------------------------------------
                    
                    IF ( h2 < hn ) THEN
                       
                       !-------------------------------------
                       ! If gap smaller than minimal 
                       ! allowed gap, set it to minimum.
                       !-------------------------------------
                       
                       IF ( h2 < hm ) THEN
                          
                          h2 = hm
                          
                       END IF
                       
                       Ft = F0 - F0*h2/hn
                       
                    END IF
                    
                 CASE (mcf_cw_repul_type_DLVO)
                    
                    !----------------------------------------
                    ! For DLVO force, it is not zero at hn
                    ! and less than 1/100 beyond 5*hn.
                    !----------------------------------------
                    
                    !----------------------------------------
                    ! If gap smaller than minimal 
                    ! allowed gap, set it to minimum.
                    !----------------------------------------
                    
                    IF ( h2 < hm ) THEN
                       
                       h2 = hm
                       
                    END IF
                    
                    Ft = F0 / hn * EXP(-h2/hn) /(1.0_MK-EXP(-h2/hn))
                    
                 END SELECT
                 
                 F(j) = F(j) - Ft
                 
                 !-------------------------------------------
                 ! Collect force on boundary
                 !-------------------------------------------
                 
                 FB(j,2*j) = Ft
                 
              END IF ! h2 < hn
              
           END IF ! bcdef(2j-1)=wall_solid
           
        END DO ! j = 1 , dim
        
        
9999    CONTINUE
        
        
        RETURN
        
      END SUBROUTINE colloid_compute_repulsion_cw
      
      
      
