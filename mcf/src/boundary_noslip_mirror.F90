      SUBROUTINE boundary_noslip_mirror(this,vw,sid_w,stat_info)
        !------------------------------------------
        !   Subroutine  : Implementing
        !                 no slip condition using
        !                 symmetry/mirror particles
        !                 created by PPM.
        !
        !  Revision     : V0.1 29.10.2009,
        !                 original version.
        !------------------------------------------
        ! Author        : Xin Bian
        ! Contact       : xin.bian@aer.mw.tum.de
        !------------------------------------------
        
        !------------------------------------------
        ! Arguments
        !------------------------------------------
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:), INTENT(INOUT)   :: vw
        INTEGER, INTENT(IN)                     :: sid_w
        
        INTEGER, INTENT(OUT)                    :: stat_info 
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        INTEGER                                 :: wall_index
        
        stat_info     = 0
        stat_info_sub = 0
        
        dim = this%num_dim
        wall_index = ABS(sid_w)
        
        vw(1:dim) = -vw(1:dim) + &
             2.0_MK*this%shear_v(1:dim,wall_index)
        
        RETURN
        
      END SUBROUTINE boundary_noslip_mirror
      
