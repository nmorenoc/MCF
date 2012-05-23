      SUBROUTINE colloid_compute_repulsion_cc(this,&
           x_ip,x_jp,sid_ip,sid_jp,F_ij,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_repulsion_cc
        !----------------------------------------------------
        !
        ! Purpose     : The gap between near contacting
        !               colloid can be too small and
        !               tend to overlap.
        !               Compute a extra repulsive force
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
        !
        ! Remarks     : V0.4 18.5.2012, shift 2nd type
        !               repulsive force down to have a zero
        !               value at 5*cut off, i.e., 5*hn.
        !
        !               V0.3 6.3.2012, change the second
        !               repulive force,assuming radius
        !               of colloid is universally one.
        !             
        !               V0.2 16.11.2010, one pair version.
        !
        !               V0.1 15.10 2010, original version,
        !               loop over all pairs
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
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:),INTENT(IN)       :: x_ip
        REAL(MK), DIMENSION(:),INTENT(IN)       :: x_jp
        INTEGER, INTENT(IN)                     :: sid_ip
        INTEGER, INTENT(IN)                     :: sid_jp
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: F_ij
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: dim,num
        REAL(MK)                        :: hn,hm, F0
        REAL(MK)                        :: aa, r, h
        REAL(MK)                        :: F
        REAL(MK), DIMENSION(3)          :: R12
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        
        dim = this%num_dim
        num = this%num_colloid
        
        F_ij(1:dim) = 0.0_MK
        F  = 0.0_MK
                
        hn = this%cc_repul_cut_off
        hm = this%cc_repul_cut_on
        F0 = this%cc_repul_F0
        
        !----------------------------------------------------
        ! Calculate the gap.
        !----------------------------------------------------

        R12(1:dim) = x_ip(1:dim) - x_jp(1:dim)
        
        r = SQRT(DOT_PRODUCT(R12(1:dim), R12(1:dim)))
        aa = this%radius(1,sid_ip) + this%radius(1,sid_jp)
        
        h = r - aa
        
        !----------------------------------------------------
        ! If gap smaller than repulsion cut off.
        !----------------------------------------------------
        
        IF ( h < 5.0_MK * hn ) THEN
           
           SELECT CASE (this%cc_repul_type)
              
           CASE (mcf_cc_repul_type_Hookean)
              !----------------------------------------------
              ! For linear spring force, it has clear zero 
              ! at hn.
              !----------------------------------------------
              
              IF ( h < hn ) THEN
                 
                 !-------------------------------------------
                 ! If gap smaller than minimal allowed gap, 
                 ! set it to minimum.
                 !-------------------------------------------
                 
                 IF ( h < hm ) THEN
                    
                    h = hm
                    
                 END IF
                 
                 F = F0 - F0*h/hn
                 
              END IF ! h < hn
              
           CASE (mcf_cc_repul_type_DLVO)
              !----------------------------------------------
              ! For DLVO force, it does not have clear zero 
              ! at hn. But at 5*hn, its value is smaller
              ! than F0/100.
              !----------------------------------------------
              
              !----------------------------------------------
              ! If gap smaller than minimal allowed gap, 
              ! set it to minimum.
              !----------------------------------------------
              
              IF ( h < hm ) THEN
                 
                 h = hm
                 
              END IF
              
              F = F0 / hn * EXP(-h/hn) /(1.0_MK-EXP(-h/hn)) - &
                   F0 / hn * EXP(-5.0_MK) /(1.0_MK-EXP(-5.0_MK))
              
           CASE DEFAULT
              
              PRINT *, __FILE__, __LINE__, &
                   "no such repulsive force!"
              stat_info = -1
              GOTO 9999
              
           END SELECT
           
           R12(1:dim) = R12(1:dim) / r        
           F_ij(1:dim) = F * R12(1:dim)
           
        END IF ! h < 5*hn
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE colloid_compute_repulsion_cc
      
      
      
