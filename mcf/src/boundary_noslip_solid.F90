      SUBROUTINE boundary_noslip_solid(this,xf,xw,vf,vw,sid_w,stat_info)
        !----------------------------------------------------
        ! Subroutine : Return the artificial velocity for a 
        !              numerical particle inside a wall,
        !              in order to get no slip condition 
        !              on the surface of the wall.
        !
        ! Remark     : Wall particles are crated by
        !              MCF using solid particles.
        !
        ! Revision   : V0.1 29.09.2009, original version
        !----------------------------------------------------
        ! Author     : Xin Bian
        ! Contact    : xin.bian@aer.mw.tum.de
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.      
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !
        ! xf    :  Position of a fulid particle.
        ! xw    :  Position of a boundary particle
        !          (inside wall).
        ! vf    :  Velocity of a fuild particle.
        ! vw    :  Artificial velocity of a 
        !          boundary particle.
        ! sid_w :  species ID of wall particle.
        !----------------------------------------------------
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xf        
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xw
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vf
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: vw
        INTEGER, INTENT(IN)                     :: sid_w
        INTEGER, INTENT(OUT)                    :: stat_info 
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: index_w
        
        
        !----------------------------------------------------
        ! Initialization.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        vw(:) = 0.0_MK
        
        index_w = ABS(sid_w)

        !----------------------------------------------------
        ! Decide which no slip type to choose.
        !----------------------------------------------------
        
        SELECT CASE( this%noslip_type )
           
        CASE (:1)
           
           CALL boundary_noslip_solid_Frozen(this,vw,index_w,stat_info_sub)
           
           IF( stat_info_sub /= 0 ) THEN
              PRINT *, "boundary_noslip : ", &
                   "Frozen boundary particle has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
        CASE (2:)
           
           CALL boundary_noslip_solid_Morris(this,xf,xw,vf,vw,index_w,&
                stat_info_sub)
           
           IF( stat_info_sub /= 0 ) THEN
              PRINT *, "boundary_noslip : ", &
                   "Morris boundary particle has problem ! "
              stat_info = -1
              GOTO 9999
           END IF
           
        END SELECT
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE boundary_noslip_solid
      
      
      SUBROUTINE boundary_noslip_solid_Frozen(this,vw,index_w,stat_info)
        !----------------------------------------------------
        !   Subroutine  : Private.
        !                 Freeze the wall boundary
        !                 particle to mimic the 
        !                 noslip condition on the
        !                 surface of the wall.
        !
        !  Revision     : V0.1 29.09.2009, 
        !                 original version.
        !----------------------------------------------------
        ! Author        : Xin Bian
        ! Contact       : xin.bian@aer.mw.tum.de
        !----------------------------------------------------

        !----------------------------------------------------
        ! Arguments
        !----------------------------------------------------
      
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: vw
        INTEGER, INTENT(IN)                     :: index_w
        INTEGER, INTENT(OUT)                    :: stat_info 
        
        
        stat_info = 0
        
        !----------------------------------------------------
        ! Assign the wall boundary particle with
        ! the speed of the wall object
        !----------------------------------------------------
        
        vw(1:this%num_dim) = this%shear_v(1:this%num_dim,index_w)
        
        RETURN
        
      END SUBROUTINE boundary_noslip_solid_Frozen
      
      
      SUBROUTINE boundary_noslip_solid_Morris(this,&
           xf,xw,vf,vw,index_w,stat_info)
        !----------------------------------------------------
        !   Subroutine  : Implementing
        !                 no slip condition from 
        !                 Morris J.P. et al. 1997.
        !
        !  Revision     : V0.1 29.09.2009,
        !                 original version.
        !----------------------------------------------------
        ! Author        : Xin Bian
        ! Contact       : xin.bian@aer.mw.tum.de
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !----------------------------------------------------
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xf
        REAL(MK), DIMENSION(:), INTENT(IN)      :: xw
        REAL(MK), DIMENSION(:), INTENT(IN)      :: vf
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: vw
        INTEGER, INTENT(IN)                     :: index_w
        INTEGER, INTENT(OUT)                    :: stat_info 
        

        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        INTEGER                                 :: index
        REAL(MK)                                :: x_sur
        REAL(MK)                                :: f_sur
        REAL(MK)                                :: w_sur
        REAL(MK)                                :: corr
        
        !----------------------------------------------------
        ! Initialization.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        dim = this%num_dim
        
        
        IF(MOD(index_w,2) == 1) THEN
           index = (index_w+1) / 2
           x_sur = this%min_phys(index)
        ELSE
           index = index_w / 2
           x_sur = this%max_phys(index)
        END IF
        
        f_sur = xf(index) - x_sur
        w_sur = xw(index) - x_sur

        !----------------------------------------------------
        ! Penetration happens
        !----------------------------------------------------
        
        IF ( f_sur * w_sur >= 0.0_MK )  THEN
           
           vw(1:dim) =  this%shear_v(1:dim,index_w)
           
        ELSE
           
           f_sur = ABS(f_sur)        
           w_sur = ABS(w_sur)
           
           !-------------------------------------------------
           ! Ratio of the distances,  used for extrapolation.
           !-------------------------------------------------
           
           corr = w_sur/f_sur
           
           !-------------------------------------------------
           ! Set the maximum ration of the extrapolation.
           !-------------------------------------------------
           
           IF (corr > mcf_wall_dist_ratio) THEN
              corr =  mcf_wall_dist_ratio
           END IF
           
           !-------------------------------------------------
           ! Extrapolated velocity for the wall 
           ! boundary particle, considering the movment 
           ! of the wall itself.
           !-------------------------------------------------
           
           vw(1:dim) = -corr * &
                ( vf(1:dim) - this%shear_v(1:dim,index_w) ) + &
                this%shear_v(1:dim,index_w)
           
        END IF
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE boundary_noslip_solid_Morris 
  
