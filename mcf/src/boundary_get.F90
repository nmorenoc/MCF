!------------------------------------------------------------
! Here are get functions/subroutines of Class Boundary
!------------------------------------------------------------

      INTEGER FUNCTION boundary_get_num_dim(this, stat_info)
        
        !---------------------------------------
        ! Return the num of dimension.
        !---------------------------------------
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_num_dim = this%num_dim
        
        RETURN
        
      END FUNCTION boundary_get_num_dim
      
      
      SUBROUTINE boundary_get_bcdef(this,d_bcdef,stat_info)
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, DIMENSION(:), POINTER  :: d_bcdef
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_bcdef)) THEN
           DEALLOCATE(d_bcdef)
        END IF
        
        ALLOCATE(d_bcdef(2*this%num_dim))
        
        d_bcdef(1:2*this%num_dim) = this%bcdef(1:2*this%num_dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_bcdef
      
      
      SUBROUTINE boundary_get_shear_rate(this,d_shear_rate,stat_info)
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_shear_rate
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_shear_rate)) THEN
           DEALLOCATE(d_shear_rate)
        END IF
        
        ALLOCATE(d_shear_rate(this%num_dim,this%num_dim))
        
        d_shear_rate(1:this%num_dim,1:this%num_dim) = &
             this%shear_rate(1:this%num_dim,1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_shear_rate

      
      SUBROUTINE boundary_get_shear_length(this,d_shear_length,stat_info)
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_shear_length
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_shear_length)) THEN
           DEALLOCATE(d_shear_length)
        END IF
        
        ALLOCATE(d_shear_length(this%num_dim,2*this%num_dim))
        
        d_shear_length(1:this%num_dim,1:2*this%num_dim) = &
             this%shear_length(1:this%num_dim,1:2*this%num_dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_shear_length

      
      SUBROUTINE boundary_get_shear_v0(this,d_shear_v0,stat_info)
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_shear_v0
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_shear_v0)) THEN
           DEALLOCATE(d_shear_v0)
        END IF
        
        ALLOCATE(d_shear_v0(this%num_dim,2*this%num_dim))
        
        d_shear_v0(1:this%num_dim,1:2*this%num_dim) = &
             this%shear_v0(1:this%num_dim,1:2*this%num_dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_shear_v0
      
      
      SUBROUTINE boundary_get_shear_v(this,d_shear_v,stat_info)
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_shear_v
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_shear_v)) THEN
           DEALLOCATE(d_shear_v)
        END IF
        
        ALLOCATE(d_shear_v(this%num_dim,2*this%num_dim))
        
        d_shear_v(1:this%num_dim,1:2*this%num_dim) = &
             this%shear_v(1:this%num_dim,1:2*this%num_dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_shear_v


      SUBROUTINE boundary_get_shear_type(this,d_shear_type,stat_info)
      
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, DIMENSION(:), POINTER  :: d_shear_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_shear_type)) THEN
           DEALLOCATE(d_shear_type)
        END IF
        
        ALLOCATE(d_shear_type(2*this%num_dim))
        
        d_shear_type(1:2*this%num_dim) = this%shear_type(1:2*this%num_dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_shear_type
      

      SUBROUTINE boundary_get_shear_freq(this,d_shear_freq,stat_info)
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:), POINTER         :: d_shear_freq
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_shear_freq)) THEN
           DEALLOCATE(d_shear_freq)
        END IF
        
        ALLOCATE(d_shear_freq(2*this%num_dim))
        
        d_shear_freq(1:2*this%num_dim) = &
             this%shear_freq(1:2*this%num_dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_shear_freq
      

      INTEGER FUNCTION boundary_get_wall_rho_type(this,stat_info)
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_wall_rho_type = this%rho_type
        
        RETURN
        
      END FUNCTION boundary_get_wall_rho_type
      

      INTEGER FUNCTION boundary_get_wall_noslip_type(this,stat_info)
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_wall_noslip_type = &
             this%noslip_type
        
        RETURN
        
      END FUNCTION   boundary_get_wall_noslip_type
      
      
      REAL(MK) FUNCTION boundary_get_dout(this,stat_info)
        
        TYPE(Boundary), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_dout = this%dout
        
        RETURN
        
      END FUNCTION boundary_get_dout
      

      SUBROUTINE boundary_get_drag(this,d_drag,stat_info)
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        
        stat_info = 0
        dim = this%num_dim
        
        IF(ASSOCIATED(d_drag)) THEN
           DEALLOCATE(d_drag)
        END IF
        
        ALLOCATE(d_drag(dim,2*dim))
        
        d_drag(1:dim,1:2*dim) = this%drag(1:dim,1:2*dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_drag
      
#ifdef __WALL_FORCE_SEPARATE
      SUBROUTINE boundary_get_drag_p(this,d_drag,stat_info)
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        
        stat_info = 0
        dim = this%num_dim
        
        IF(ASSOCIATED(d_drag)) THEN
           DEALLOCATE(d_drag)
        END IF
        
        ALLOCATE(d_drag(dim,2*dim))
        
        d_drag(1:dim,1:2*dim) = this%drag_p(1:dim,1:2*dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_drag_p
      
      
      SUBROUTINE boundary_get_drag_v(this,d_drag,stat_info)
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        
        stat_info = 0
        dim = this%num_dim
        
        IF(ASSOCIATED(d_drag)) THEN
           DEALLOCATE(d_drag)
        END IF
        
        ALLOCATE(d_drag(dim,2*dim))
        
        d_drag(1:dim,1:2*dim) = this%drag_v(1:dim,1:2*dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_drag_v
      
      
      SUBROUTINE boundary_get_drag_r(this,d_drag,stat_info)
        
        TYPE(Boundary), INTENT(IN)              :: this
        REAL(MK), DIMENSION(:,:), POINTER       :: d_drag
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: dim
        
        stat_info = 0
        dim = this%num_dim
        
        IF(ASSOCIATED(d_drag)) THEN
           DEALLOCATE(d_drag)
        END IF
        
        ALLOCATE(d_drag(dim,2*dim))
        
        d_drag(1:dim,1:2*dim) = this%drag_r(1:dim,1:2*dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_drag_r

#endif
      
      SUBROUTINE boundary_get_min_phys(this,d_min_phys,stat_info)
        
        TYPE(Boundary), INTENT(IN)      :: this
        REAL(MK), DIMENSION(:), POINTER :: d_min_phys
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_min_phys)) THEN
           DEALLOCATE(d_min_phys)
        END IF
        
        ALLOCATE(d_min_phys(this%num_dim))
        
        d_min_phys(1:this%num_dim) = this%min_phys(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_min_phys
      
      
      SUBROUTINE boundary_get_max_phys(this,d_max_phys,stat_info)
        
        TYPE(Boundary), INTENT(IN)      :: this
        REAL(MK), DIMENSION(:), POINTER :: d_max_phys
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_max_phys)) THEN
           DEALLOCATE(d_max_phys)
        END IF
        
        ALLOCATE(d_max_phys(this%num_dim))
        
        d_max_phys(1:this%num_dim) = this%max_phys(1:this%num_dim)
        
        RETURN
        
      END SUBROUTINE  boundary_get_max_phys
      

      INTEGER FUNCTION boundary_get_num_peri(this, stat_info)
        
        !---------------------------------------
        ! Return the num of periodic boundaires
        !---------------------------------------
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_num_peri = this%num_peri
        
        RETURN
        
      END FUNCTION boundary_get_num_peri
      

      INTEGER FUNCTION boundary_get_num_sym(this, stat_info)
        
        !---------------------------------------
        ! Return the num of symmetry boundaires
        !---------------------------------------
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_num_sym = this%num_sym
        
        RETURN
        
      END FUNCTION boundary_get_num_sym


      INTEGER FUNCTION boundary_get_num_wall_sym(this, stat_info)
        
        !---------------------------------------
        ! Return the num of wall boundaires,
        ! handeled by PPM.
        !---------------------------------------
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_num_wall_sym = this%num_wall_sym
        
        RETURN
        
      END FUNCTION boundary_get_num_wall_sym


      INTEGER FUNCTION boundary_get_num_wall_solid(this, stat_info)
        
        !---------------------------------------
        ! Return the num of wall boundaires,
        ! handeled by MCF.
        !---------------------------------------
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_num_wall_solid = this%num_wall_solid
        
        RETURN
        
      END FUNCTION boundary_get_num_wall_solid
      
      
      INTEGER FUNCTION boundary_get_num_wall(this, stat_info)
        
        !---------------------------------------
        ! Return the num of wall boundaires,
        ! handeled either by PPM, or MCF.
        !---------------------------------------
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_num_wall = &
             this%num_wall_sym + this%num_wall_solid
        
        RETURN
        
      END FUNCTION boundary_get_num_wall
      
      
      INTEGER FUNCTION boundary_get_num_le(this, stat_info)
        
        !---------------------------------------
        ! Return the num of Lees-Edwards boundaires.
        !---------------------------------------
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_num_le = this%num_le
        
        RETURN
        
      END FUNCTION boundary_get_num_le
      
      
      INTEGER FUNCTION boundary_get_num_shear(this, stat_info)
        
        !---------------------------------------
        ! Return the num of shear boundaires.
        !---------------------------------------
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_num_shear = &
             this%num_shear
        
        RETURN
        
      END FUNCTION boundary_get_num_shear


      INTEGER FUNCTION boundary_get_num_part_wall_solid(this, stat_info)
        
        !---------------------------------------
        ! Return the num of wall bounday particles
        ! created by MCF.
        !---------------------------------------
        
        TYPE(Boundary), INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        boundary_get_num_part_wall_solid = this%num_part_wall_solid
        
        RETURN
        
      END FUNCTION boundary_get_num_part_wall_solid
      
    
      SUBROUTINE boundary_get_wall_particle_v(this,num_dim,vw,sid,stat_info)
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
        INTEGER, INTENT(IN)                     :: num_dim
        REAL(MK), DIMENSION(:), INTENT(OUT)     :: vw
        INTEGER, INTENT(IN)                     :: sid
        INTEGER, INTENT(OUT)                    :: stat_info 
        
        INTEGER                                 :: index_w
        
        stat_info = 0
        
        index_w = ABS(sid)

        IF ( num_dim /= this%num_dim ) THEN
           
           PRINT *, "boundary_get_wall_particle_v : ", &
                "Dimension doesn't match !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Assign the wall boundary particle with
        ! the speed of the wall object
        !----------------------------------------------------
        
        vw(1:num_dim) = this%shear_v(1:num_dim,index_w)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE boundary_get_wall_particle_v
      
