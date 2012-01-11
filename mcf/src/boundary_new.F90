!-------------------------------------------------
! Constructors of Class Boundary
!-------------------------------------------------
      SUBROUTINE boundary_init_default(this,stat_info)
        !------------------------------------------
        ! Default constructor, set to 2D.
        !------------------------------------------
        
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%num_dim = 2
        
        NULLIFY(this%bcdef)
        ALLOCATE(this%bcdef(4)) 
        this%bcdef(1:4) = ppm_param_bcdef_periodic

        NULLIFY(this%shear_rate)
        ALLOCATE(this%shear_rate(2,2))
        this%shear_rate(1:2,1:2) = 0.0_MK
        
        NULLIFY(this%shear_length)
        ALLOCATE(this%shear_length(2,4))
        this%shear_length(1:2,1:4) = 0.0_MK
        
        NULLIFY(this%shear_type)
        ALLOCATE(this%shear_type(4))
        this%shear_type(1:4) = 0
        
        NULLIFY(this%shear_v0)
        ALLOCATE(this%shear_v0(2,4))
        this%shear_v0(1:2,1:4) = 0.0_MK
        
        NULLIFY(this%shear_v)
        ALLOCATE(this%shear_v(2,4))
        this%shear_v(1:2,1:4) = 0.0_MK
        
        NULLIFY(this%shear_freq)
        ALLOCATE(this%shear_freq(4))
        this%shear_freq(1:4) = 0.0_MK
        
        this%rho_type = 0

        this%noslip_type = 1
        
        this%dout = 0.0_MK
        
        NULLIFY(this%drag)
        ALLOCATE(this%drag(2,4))
        this%drag(:,:) = 0.0_MK
#ifdef __FORCE_SEPARATE
        NULLIFY(this%drag_p)
        ALLOCATE(this%drag_p(2,4))
        this%drag_p(:,:) = 0.0_MK
        NULLIFY(this%drag_v)
        ALLOCATE(this%drag_v(2,4))
        this%drag_v(:,:) = 0.0_MK
        NULLIFY(this%drag_r)
        ALLOCATE(this%drag_r(2,4))
        this%drag_r(:,:) = 0.0_MK
#endif
        
        this%min_phys(:)     = 0.0_MK
        this%max_phys(1:2)   = 1.0e-3_MK
        this%min_phys_t(:)   = 0.0_MK        
        this%max_phys_t(1:2) = 1.0e-3_MK
        
        this%num_peri       = 0
        this%num_sym        = 0
        this%num_wall_sym   = 0
        this%num_wall_solid = 0
        this%num_osci       = 0
        this%num_le         = 0
        this%num_shear      = 0
        
        this%num_part_wall_solid = 0
        
        RETURN
        
      END SUBROUTINE boundary_init_default

      
      SUBROUTINE boundary_init(this,d_num_dim,stat_info)

        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER, INTENT(IN)             :: d_num_dim
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%num_dim = d_num_dim
        
        NULLIFY(this%bcdef)
        ALLOCATE(this%bcdef(2*d_num_dim)) 
        this%bcdef(1:2*d_num_dim) = 0
    
        NULLIFY(this%shear_rate)
        ALLOCATE(this%shear_rate(d_num_dim,d_num_dim))
        this%shear_rate(1:d_num_dim,1:d_num_dim) = 0.0_MK
        
        NULLIFY(this%shear_length)
        ALLOCATE(this%shear_length(d_num_dim,2*d_num_dim))
        this%shear_length(1:d_num_dim,1:2*d_num_dim) = 0.0_MK
       
        NULLIFY(this%shear_type)
        ALLOCATE(this%shear_type(2*d_num_dim))
        this%shear_type(1:2*d_num_dim) = 0
        
        NULLIFY(this%shear_v0)
        ALLOCATE(this%shear_v0(d_num_dim,2*d_num_dim))
        this%shear_v0(1:d_num_dim,1:2*d_num_dim) = 0.0_MK
        
        NULLIFY(this%shear_v)
        ALLOCATE(this%shear_v(d_num_dim,2*d_num_dim))
        this%shear_v(1:d_num_dim,1:2*d_num_dim) = 0.0_MK
        
        NULLIFY(this%shear_freq)
        ALLOCATE(this%shear_freq(2*d_num_dim))
        this%shear_freq(1:2*d_num_dim) = 0.0_MK
    
        this%rho_type = 0
        
        this%noslip_type = 1
        
        this%dout = 0.0_MK
        
        NULLIFY(this%drag)
        ALLOCATE(this%drag(d_num_dim,2*d_num_dim))
        this%drag(:,:) = 0.0_MK
#ifdef __FORCE_SEPARATE
        NULLIFY(this%drag_p)
        ALLOCATE(this%drag_p(d_num_dim,2*d_num_dim))
        this%drag_p(:,:) = 0.0_MK
        NULLIFY(this%drag_v)
        ALLOCATE(this%drag_v(d_num_dim,2*d_num_dim))
        this%drag_v(:,:) = 0.0_MK
        NULLIFY(this%drag_r)
        ALLOCATE(this%drag_r(d_num_dim,2*d_num_dim))
        this%drag_r(:,:) = 0.0_MK
#endif

        this%min_phys(:) = 0.0_MK
        this%max_phys(:) = 0.0_MK
        this%min_phys_t(:) = 0.0_MK
        this%max_phys_t(:) = 0.0_MK
        
        this%num_peri       = 0 
        this%num_sym        = 0
        this%num_wall_sym   = 0
        this%num_wall_solid = 0
        this%num_osci       = 0
        this%num_le         = 0 
        this%num_shear      = 0 
                
        this%num_part_wall_solid = 0
        
        RETURN
        
      END SUBROUTINE boundary_init
      
      
      SUBROUTINE boundary_display_parameters(this,stat_info)
        
        !----------------------------------------------------
        ! Subroutine : boundary_display_parameters
        !
        !
        !----------------------------------------------------

        TYPE(Boundary),INTENT(IN)       :: this
        INTEGER,INTENT(OUT)             :: stat_info
        
        INTEGER                         :: dim,j
        INTEGER                         :: num_wall
        INTEGER                         :: num_shear
        INTEGER                         :: stat_info_sub

        stat_info = 0
        stat_info_sub = 0
        
        dim =this%num_dim
        num_wall  = boundary_get_num_wall(this,stat_info_sub)
        num_shear = boundary_get_num_shear(this,stat_info_sub)
        
        PRINT *, '---***************Start***************---'
        PRINT *, '     Boundary parameters'
        PRINT *, '---***********************************---'
        
        PRINT *, "num_dim          : ", dim
        PRINT *, "bcdef            : "
        PRINT *,  this%bcdef(1:2*dim)
        
        IF ( num_shear > 0 ) THEN
           
           PRINT *, "shear type       : "
           PRINT *, this%shear_type(1:2*dim)
           PRINT *, "shear rate       : "
           IF (dim==2) THEN
              PRINT *, this%shear_rate(2,1)
              PRINT *, this%shear_rate(1,2)
           ELSE IF (dim ==3 ) THEN
              PRINT *, this%shear_rate(2:3,1)
              PRINT *, this%shear_rate(1,2),this%shear_rate(3,2)
              PRINT *, this%shear_rate(1,3),this%shear_rate(2,3)
           END IF
           PRINT *, "shear velocity   : " 
           
           IF (dim==2) THEN
              PRINT *, this%shear_v0(2,1:2)
              PRINT *, this%shear_v0(1,3:4)
           ELSE IF (dim ==3 ) THEN
              PRINT *, this%shear_v0(2:3,1)
              PRINT *, this%shear_v0(2:3,2)
              PRINT *, this%shear_v0(1,3),this%shear_v0(3,3)
              PRINT *, this%shear_v0(1,4),this%shear_v0(3,4)
              PRINT *, this%shear_v0(1:2,5)
              PRINT *, this%shear_v0(1:2,6)
           END IF
           
           PRINT *, "shear frequency  : "
           DO j=1,dim
              PRINT *, this%shear_freq(j*2-1:j*2)
           END DO
           
           IF ( this%num_wall_solid > 0 ) THEN
              
              PRINT *, "wall rho type    : ", this%rho_type
              
           END IF

           IF ( num_wall > 0 ) THEN
              PRINT *, "no slip          : "
              PRINT *, this%noslip_type
           END IF
           
        END IF
        
        
        PRINT *, '---****************End****************---'
        
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE boundary_display_parameters
      
      
