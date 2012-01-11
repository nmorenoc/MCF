      SUBROUTINE colloid_check_penetration(this, &
           xf_in,vf_in,sid_c,vf_out,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_check_penetration
        !----------------------------------------------------
        ! Purpose     : 
        !
        ! Revision    : V0.1 28.04.2010, orignal version.
        !               Check penetration, and use
        !               bounce back rule for fluid. 
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
        ! xf    : position of a fluid particle.
        ! vf_in : velocity of a fluid particle.
        ! sid_c : species ID of a colloid boundary particle.
        !
        ! Output
        !
        ! vf_out : velocity after bounce back for fluid particle.
        ! stat_info : status of the routine.
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(IN)               :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xf_in
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vf_in
        INTEGER, INTENT(IN)                     :: sid_c
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: vf_out
        INTEGER, INTENT(OUT)                    :: stat_info 
        
        !----------------------------------------------------
        ! Local variables start here :
        ! 
        ! i     : index.
        ! dim   : number of dimension.
        ! xcoll : position of the center of the colloid.
        ! vcoll : velocity of the center of the colloid.
        ! num_le: number of Lees-Edwards bondaries.
        ! shear_length :
        ! shear_v :
        ! length  : box length.
        ! r_xf    : relative position from fluid particle 
        !           to the center.
        ! r_xc    : relative position from colloid boundary
        !           particle to center.
        ! r_vc    : relative velocity of colloid boundary 
        !           particle to the center.
        ! d_fcoll : distance from fluid particle to the center.
        ! d_ccoll : distance from the colloid boundary particle
        !           to the center.
        ! nvector : unit vector pointing from center to
        !           fluid particle.
        ! d_cn    : distance mapping d_ccoll to nvector.
        ! d_fs    : distance from fluid to tangent surface.
        ! d_cs    : distance from fluid to tangent surface.
        ! corr    : d_cs/d_fs.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: i
        INTEGER                                 :: dim
        INTEGER                                 :: num_le
        REAL(MK),DIMENSION(:,:),POINTER         :: shear_length
        REAL(MK),DIMENSION(:,:),POINTER         :: shear_v
        REAL(MK), DIMENSION(3)                  :: length
        REAL(MK), DIMENSION(3)                  :: xcoll
        REAL(MK), DIMENSION(3)                  :: vcoll
        REAL(MK), DIMENSION(3)                  :: r_xf
        REAL(MK)                                :: d_fcoll
        REAL(MK)                                :: corr
          
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        dim = this%num_dim
        NULLIFY(shear_length)
        NULLIFY(shear_v)
        
        !----------------------------------------------------
        ! Get shear length, velocity for Lees-Edwards boundary.
        !----------------------------------------------------
        
        num_le = &
             boundary_get_num_le(this%boundary,stat_info_sub)
        
        IF ( num_le > 0 ) THEN
           
           CALL boundary_get_shear_length(this%boundary,&
                shear_length,stat_info_sub)
           CALL boundary_get_shear_v(this%boundary,&
                shear_v,stat_info_sub)

        END IF
        
        !----------------------------------------------------
        ! Calculate the length of computational box.
        !----------------------------------------------------
        
        length(1:dim) = &
             this%max_phys(1:dim) - this%min_phys(1:dim)
        
        !----------------------------------------------------
        ! Get the colloid center's position and velocity.
        !----------------------------------------------------
        
        xcoll(1:dim) = this%x(1:dim,sid_c)
        vcoll(1:dim) = this%v(1:dim,sid_c)
        
        !----------------------------------------------------
        ! Get relative position of fluid and colloid boundary
        ! particle to the center.
        !----------------------------------------------------
        
        r_xf(1:dim) = xf_in(1:dim) - xcoll(1:dim)
        
        !----------------------------------------------------
        ! In context of periodic or Lees-Edwards boundary, 
        ! we have to consider the colloid boundary particle 
        ! which may be inside the image of colloid, instead 
        ! of colloid itself.
        !----------------------------------------------------
        
        Do i = 1, dim
           
           !-------------------------------------------------
           ! Check if the colloid boundary particle is
           ! around the colloid center.
           ! If yes, use the center coordinate of colloid.
           ! If no, translate the center coodinate
           ! of the colloid to get the image center.
           !-------------------------------------------------
           
           IF( r_xc(i) > this%ra(sid_c) + this%dout ) THEN
              
              SELECT CASE ( this%bcdef(2*i-1) ) 
                 
              CASE ( ppm_param_bcdef_periodic )
                 
                 xcoll(i) = xcoll(i) + length(i)
                 
              CASE ( ppm_param_bcdef_LE )
                 
                 DO k = 1, dim
                    
                    IF ( k == i ) THEN
                       
                       xcoll(k) = xcoll(k) + length(k)
                       
                    ELSE
                       
                       xcoll(k) = &
                            MODULO(xcoll(k) - &
                            shear_length(k,2*i-1), &
                            length(k) )
                       
                       vcoll(k) = vcoll(k) + &
                            ( shear_v(k,2*i) - shear_v(k,2*i-1) )
                       
                    END IF
                    
                 END DO ! k
                 
              END SELECT ! bcdef
              
           ELSE IF( -r_xc(i) > this%ra(sid_c) + this%dout ) THEN
              
              SELECT CASE ( this%bcdef(2*i) )
                 
              CASE ( ppm_param_bcdef_periodic )
                 
                 xcoll(i) = xcoll(i) - length(i)
                 
              CASE ( ppm_param_bcdef_LE )
                 
                 DO k = 1, dim
                    
                    IF ( k == i ) THEN
                       
                       xcoll(k) = xcoll(k) - length(k)
                       
                    ELSE
                       
                       xcoll(k) = &
                            MODULO(xcoll(k) - &
                            shear_length(k,2*i), &
                            length(k) )
                       
                       vcoll(k) = vcoll(k) + &
                            ( shear_v(k,2*i-1) - shear_v(k,2*i) )
                       
                    END IF
                    
                 END DO ! k
                 
              END SELECT ! bcdef
              
           END IF ! r_xc(i) > radius + dout
           
        END DO ! i
        
        !----------------------------------------------------
        ! Get relative position of fluid particle 
        ! to the center again.
        !----------------------------------------------------
        
        r_xf(1:dim) = xf_in(1:dim) - xcoll(1:dim)
        
        !----------------------------------------------------
        ! Distance between a fluid particle and the center.
        !----------------------------------------------------
        
        d_fcoll = 0.0_MK
        DO i = 1, dim
           d_fcoll = d_fcoll + r_xf(i)**2
        END DO
        d_fcoll = SQRT(d_fcoll)
        
        !----------------------------------------------------
        ! Check if penetration happens, i.e., a fluid 
        ! particle goes inside a colloid object. If yes, 
        ! simply assign the same velocity as the colloid.
        !----------------------------------------------------
        
        IF( d_fcoll <= this%ra(sid_c) ) THEN
           
           vf_out(1:dim) = - vf_in(1:dim)
           
        END IF
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Release dynamic memory.
        !----------------------------------------------------
        
        IF(ASSOCIATED(shear_length)) THEN
           DEALLOCATE(shear_length) 
        END IF
        
        IF(ASSOCIATED(shear_v)) THEN
           DEALLOCATE(shear_v) 
        END IF
        
        RETURN
        
      END SUBROUTINE colloid_check_penetration
