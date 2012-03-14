      SUBROUTINE colloid_compute_lubrication_cc(this,&
           x_ip,x_jp,v_ip,v_jp,omega_ip,omega_jp,&
           sid_ip,sid_jp,F_i,F_j,T_i,T_j,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_compute_lubrication_cc
        !----------------------------------------------------
        !
        ! Purpose     : When the resolution is finite, fluid 
        !               in the gap between near contacting
        !               colloid can not be resolved. 
        !               We use lubrication theory to
        !               compute the lubrication force.
        !
        ! Routines    :
        !
        ! References  :
        !               If lubrication correction is made
        !               to colloid-colloid interaction,
        !               2D and 3D implementations are
        !               performed.
        !               (2D : Kromkamp et al. Chem.Eng.Sci.2006)
        !               (3D : Nguyen and Ladd, Phys.Rev.E 2002)
        !
        !
        ! Remarks     : V0.3 27.1.2011, tangential force and
        !               torque are added for equal-sized cylinders.
        !
        !               V0.2 16.11.2010, one pair version.
        !              
        !               V0.1 27.08 2010, original version.
        !               loop over all pairs.
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
        REAL(MK), DIMENSION(:),INTENT(IN)       :: v_ip
        REAL(MK), DIMENSION(:),INTENT(IN)       :: v_jp
        REAL(MK), DIMENSION(1:3),INTENT(IN)     :: omega_ip
        REAL(MK), DIMENSION(1:3),INTENT(IN)     :: omega_jp      
        INTEGER, INTENT(IN)                     :: sid_ip
        INTEGER, INTENT(IN)                     :: sid_jp
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: F_i
        REAL(MK), DIMENSION(:),INTENT(OUT)      :: F_j
        REAL(MK), DIMENSION(1:3),INTENT(OUT)    :: T_i
        REAL(MK), DIMENSION(1:3),INTENT(OUT)    :: T_j
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        ! 
        ! F_n : normal direction force.
        ! F_t : tangential direction force.
        ! T   : torque
        ! e12 perpendicular to g12
        !----------------------------------------------------
        
        INTEGER                         :: dim
        REAL(MK)                        :: hn,hm
        REAL(MK)                        :: a, a1,a2,aa,r,h
        REAL(MK), DIMENSION(3)          :: u12,r12,e12,g12
        
        REAL(MK)                        :: F0,F1
        REAL(MK)                        :: F_n
        REAL(MK)                        :: F_t1, F_t2
        REAL(MK)                        :: T1, T2

        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        
        dim = this%num_dim
        hn  = this%cc_lub_cut_off
        hm  = this%cc_lub_cut_on
        
        F_i(1:dim) = 0.0_MK
        F_j(1:dim) = 0.0_MK
        
        T_i(1:3) = 0.0_MK
        T_j(1:3) = 0.0_MK
        
        !----------------------------------------------------
        ! Calculate the gap.
        !----------------------------------------------------
        
        r12(1:dim) = x_ip(1:dim) - x_jp(1:dim)
        r  = SQRT(DOT_PRODUCT(r12(1:dim), r12(1:dim)))
        a1 = this%radius(1,sid_ip)
        a2 = this%radius(1,sid_jp)
        a  = a1
        aa = a1 + a2
        h = r - aa
        
        u12(1:dim) = v_ip(1:dim) - v_jp(1:dim)
        e12(1:dim) = r12(1:dim) / r
        
        !----------------------------------------------------
        ! Calculate normal force.
        !----------------------------------------------------
        
        IF ( h < hn ) THEN
            
           !-------------------------------------------------
           ! If gap smaller than minimal allowed gap, 
           ! set it to minimum.
           !-------------------------------------------------
           
           IF ( h < hm ) THEN
              
              h = hm
              
           END IF
           
           F_n=0.0_MK
           
           IF ( dim == 2  ) THEN
              
              F0  = 3.0_MK*mcf_pi*SQRT(2.0_MK)/4.0_MK
              F1  = 231.0_MK*mcf_pi*SQRT(2.0_MK) / 80.0_MK
              
              F_n = -0.5_MK * this%eta * &
                   DOT_PRODUCT(u12(1:dim), e12(1:dim)) * &
                   ( (aa/h)**1.5_MK  * (F0 + h*F1/aa)  - &
                   (aa/hn)**1.5_MK * (F0 + hn*F1/aa) )
              
              
           ELSE IF (dim == 3 ) THEN
              
              F_n = -6.0_MK * mcf_pi * this%eta * &
                   DOT_PRODUCT(u12(1:dim), e12(1:dim)) * &
                   (a1*a2/aa)**2 * &
                   ( 1.0_MK/h - 1.0_MK/hn )
              
           END IF ! dim
           
           F_i(1:dim) = F_i(1:dim) + F_n * e12(1:dim)
           F_j(1:dim) = F_j(1:dim) - F_n * e12(1:dim)
           
        END IF ! h < hn

        !----------------------------------------------------
        ! Calculate tangential force.
        !----------------------------------------------------

        IF ( h < hn ) THEN
           
           !-------------------------------------------------
           ! If gap smaller than minimal allowed gap, 
           ! set it to minimum.
           !-------------------------------------------------
           
           IF ( h < hm ) THEN
              
              h = hm
              
           END IF
           
           F_t1=0.0_MK
           F_t2=0.0_MK
           
           IF ( dim == 2  ) THEN
              
              g12(1) = e12(2)
              g12(2) = -e12(1)
              
#if 0
              F_t1 = -this%eta * &
                   ( (17.0_MK * DOT_PRODUCT(v_ip(1:dim), g12(1:dim)) + &
                   DOT_PRODUCT(v_jp(1:dim), g12(1:dim) ) )/8.0_MK  * &
                   ( mcf_pi*SQRT(a/h) - mcf_pi*SQRT(a/hn) ) + &
                   a * (5.0_MK*omega_ip(3)-omega_jp(3))/2.0_MK * &
                   ( mcf_pi*SQRT(a/h) - mcf_pi*SQRT(a/hn) ) ) 
              
              F_t2 = -this%eta * &
                   ( ( DOT_PRODUCT(v_ip(1:dim), g12(1:dim)) + &
                   17.0_MK*DOT_PRODUCT(v_jp(1:dim), g12(1:dim) ) )/8.0_MK  * &
                   ( mcf_pi*SQRT(a/h) - mcf_pi*SQRT(a/hn) ) + &
                   a * (omega_ip(3)-5.0_MK*omega_jp(3))/2.0_MK * &
                   ( mcf_pi*SQRT(a/h)- mcf_pi*SQRT(a/hn) ) ) 
#endif

#if 0
              !1st and 2nd order term
              F_t1 = -this%eta * &
                   ( (17.0_MK * DOT_PRODUCT(v_ip(1:dim), g12(1:dim)) + &
                   DOT_PRODUCT(v_jp(1:dim), g12(1:dim) ) )/8.0_MK  * &
                   ( mcf_pi*SQRT(a/h)+mcf_pi*SQRT(h/a) - &
                   mcf_pi*SQRT(a/hn)-mcf_pi*SQRT(hn/a) ) + &
                   a * (5.0_MK*omega_ip(3)-omega_jp(3))/2.0_MK * &
                   ( mcf_pi*SQRT(a/h)-mcf_pi*SQRT(h/a) - &
                   mcf_pi*SQRT(a/hn)+mcf_pi*SQRT(hn/a) ) ) 
              
              F_t2 = -this%eta * &
                   ( ( DOT_PRODUCT(v_ip(1:dim), g12(1:dim)) + &
                   17.0_MK*DOT_PRODUCT(v_jp(1:dim), g12(1:dim) ) )/8.0_MK  * &
                   ( mcf_pi*SQRT(a/h)+mcf_pi*SQRT(h/a) - &
                   mcf_pi*SQRT(a/hn)-mcf_pi*SQRT(hn/a) ) + &
                   a * (omega_ip(3)-5.0_MK*omega_jp(3))/2.0_MK * &
                   ( mcf_pi*SQRT(a/h)-mcf_pi*SQRT(h/a) - &
                   mcf_pi*SQRT(a/hn)+mcf_pi*SQRT(hn/a) ) ) 
     
#endif         
           END IF ! dim
           
           F_i(1:dim) = F_i(1:dim) + F_t1 * g12(1:dim)
           F_j(1:dim) = F_j(1:dim) + F_t2 * g12(1:dim)
           
        END IF ! h < hn

        !----------------------------------------------------
        ! Calculate torque.
        !----------------------------------------------------

        IF ( h < hn ) THEN
           
           !-------------------------------------------------
           ! If gap smaller than minimal allowed gap, 
           ! set it to minimum.
           !-------------------------------------------------
           
           IF ( h < hm ) THEN
              
              h = hm
              
           END IF
           
           T1=0.0_MK
           T2=0.0_MK
        
           IF ( dim == 2  ) THEN
              
              g12(1) = e12(2)
              g12(2) = -e12(1)
              
#if 0              
              T1 = -this%eta * a *&
                   ( (5.0_MK * DOT_PRODUCT(v_ip(1:dim), g12(1:dim)) + &
                   DOT_PRODUCT(v_jp(1:dim), g12(1:dim) ) )/2.0_MK * &
                   ( mcf_pi*SQRT(a/h) - mcf_pi*SQRT(a/hn) )  + &
                   2.0_MK * a * (2.0_MK*omega_ip(3)-omega_jp(3))  * &
                   ( mcf_pi*SQRT(a/h) - mcf_pi*SQRT(a/hn) ) )
              
              T2 = this%eta * a *&
                   ( ( DOT_PRODUCT(v_ip(1:dim), g12(1:dim)) + &
                   5.0_MK *DOT_PRODUCT(v_jp(1:dim), g12(1:dim) ) )/2.0_MK * &
                   ( mcf_pi*SQRT(a/h) - mcf_pi*SQRT(a/hn) ) + &
                   2.0_MK * a * (omega_ip(3)-2.0_MK*omega_jp(3))  * &
                   ( mcf_pi*SQRT(a/h) - mcf_pi*SQRT(a/hn) ) )
#endif
              
#if 0
              !1st and 2nd order term              
              T1 = -this%eta * a *&
                   ( (5.0_MK * DOT_PRODUCT(v_ip(1:dim), g12(1:dim)) + &
                   DOT_PRODUCT(v_jp(1:dim), g12(1:dim) ) )/2.0_MK * &
                   ( mcf_pi* SQRT(a/h)-mcf_pi*SQRT(h/a) - &
                   mcf_pi*SQRT(a/hn)+mcf_pi*SQRT(hn/a) ) + &
                   2.0_MK * a * (2.0_MK*omega_ip(3)-omega_jp(3))  * &
                   ( mcf_pi* SQRT(a/h)+1.5_MK*mcf_pi*SQRT(h/a) - &
                   mcf_pi*SQRT(a/hn)-1.5_MK*mcf_pi*SQRT(hn/a) ) )
              
              T2 = this%eta * a *&
                   ( ( DOT_PRODUCT(v_ip(1:dim), g12(1:dim)) + &
                   5.0_MK *DOT_PRODUCT(v_jp(1:dim), g12(1:dim) ) )/2.0_MK * &
                   ( mcf_pi* SQRT(a/h)-mcf_pi*SQRT(h/a) - &
                   mcf_pi*SQRT(a/hn)+mcf_pi*SQRT(hn/a) ) + &
                   2.0_MK * a * (omega_ip(3)-2.0_MK*omega_jp(3))  * &
                   ( mcf_pi* SQRT(a/h)+1.5_MK*mcf_pi*SQRT(h/a) - &
                   mcf_pi*SQRT(a/hn)-1.5_MK*mcf_pi*SQRT(hn/a) ) )
#endif     
         
           END IF ! dim           
           
           T_i(3) = T_i(3) + T1
           T_j(3) = T_j(3) + T2

           
        END IF ! h < hn
        
        RETURN          
        
      END SUBROUTINE colloid_compute_lubrication_cc
      
      
      
