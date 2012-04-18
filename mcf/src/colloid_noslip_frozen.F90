      SUBROUTINE colloid_noslip_frozen(this,xc,vc,sid_c,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_noslip_frozen
        !----------------------------------------------------
        !
        ! Purpose     : Private.
        !               Freeze the colloid boundary particle
        !               to mimic the noslip condition on the 
        !               surface of the colloid
        !
        ! Revision    : V0.2 27.11.2009, including
        !               Lees-Edwards boundary.
        !
        !               V0.1 01.03.2009, original version.
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
        ! Arguments :
        !
        ! Input 
        !
        ! this  : object of a colloid.
        ! xc    : position of a colloid boundary particle.
        ! sid_c : species ID of a colloid boundary particle.
        !
        ! Output
        !
        ! vc    : extrapolated velocity for the colloid
        !         boundary particle.
        ! stat_info : status of the routine.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xc
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: vc
        INTEGER, INTENT(IN)                     :: sid_c            
        INTEGER, INTENT(OUT)                    :: stat_info 
        
        !----------------------------------------------------
        ! Local variables.
        !
        ! i     : index.
        ! dim   : number of dimension.
        ! xcoll : position of the center of the colloid.
        ! vcoll : velocity of the center of the colloid.
        ! num_le: number of Lees-Edwards bondaries.
        ! shear_length :
        ! shear_v :
        ! length  : box length.
        ! r_xc    : relative position from colloid boundary
        !           particle to center.
        ! r_vc    : relative velocity of colloid boundary 
        !           particle to the center.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        REAL(MK), DIMENSION(3)                  :: xcoll
        REAL(MK), DIMENSION(3)                  :: vcoll
        REAL(MK), DIMENSION(3)                  :: r_xc
        REAL(MK), DIMENSION(3)                  :: r_vc
        
          
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        dim  = this%num_dim
        r_xc(:) = 0.0_MK
        
        !----------------------------------------------------        
        ! Get the nearest image colloid center's to the
        ! boundary particle.
        !----------------------------------------------------
        
        CALL colloid_in_nearest_image(this,xc(1:dim),sid_c, &
             xcoll(1:dim),r_xc(1:dim),vcoll(1:dim),stat_info_sub)

        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, "colloid_noslip_Frozen: ",&
                "colloid in_nearst_image failed !"
           stat_info = -1
           GOTO 9999
        END IF
        !----------------------------------------------------
        ! Assign the colloid boundary particle with
        ! the speed of the colloid center.
        !----------------------------------------------------
        
        vc(1:dim) =  vcoll(1:dim) 
        
        !----------------------------------------------------
        ! Add up equavilent translational velocity from
        ! rotation
        !----------------------------------------------------
        
        r_vc(:) = 0.0_MK
        
        IF ( this%rotate ) THEN
           
           CALL tool_cross_product(this%tool,&
                this%omega(1:3,sid_c,1), r_xc(1:3),&
                r_vc(1:3),stat_info_sub)
           
           vc(1:dim) =  vc(1:dim) + r_vc(1:dim)
           
        END IF
        
        
9999    CONTINUE
        
        
        RETURN
        
      END SUBROUTINE colloid_noslip_frozen
      
