!------------------------------------------------------------
! Here are set functions/subroutines of Class Boundary
!------------------------------------------------------------

      SUBROUTINE boundary_set_num_dim(this,d_num_dim, stat_info)
        !-----------------------------------------
        ! Set the number of dimension of boundary.
        !-----------------------------------------
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER, INTENT(IN)             :: d_num_dim
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        !---------------------------------------
        ! Only 2D, 3D are supported
        !---------------------------------------
        
        IF(d_num_dim < 2 .OR. d_num_dim > 3 ) THEN
           PRINT *, "boundary_set_num_dim : ", &
                "Dimension is not supported !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !---------------------------------------
        ! If reset the dimension, 
        ! the memory has to be reallocated
        !---------------------------------------
        
        IF( d_num_dim /= this%num_dim ) THEN
           
           this%num_dim = d_num_dim
           
           IF (ASSOCIATED(this%bcdef)) THEN
              DEALLOCATE(this%bcdef)
           END IF
           ALLOCATE(this%bcdef(2*d_num_dim)) 
           this%bcdef(:) = 0

           IF (ASSOCIATED(this%shear_rate)) THEN
              DEALLOCATE(this%shear_rate)
           END IF
           ALLOCATE(this%shear_rate(d_num_dim,d_num_dim))
           this%shear_rate(:,:) = 0.0_MK

           IF (ASSOCIATED(this%shear_length)) THEN
              DEALLOCATE(this%shear_length)
           END IF
           ALLOCATE(this%shear_length(d_num_dim,2*d_num_dim))
           this%shear_length(:,:) = 0.0_MK
           
           IF (ASSOCIATED(this%shear_type)) THEN
              DEALLOCATE(this%shear_type)
           END IF
           ALLOCATE(this%shear_type(2*d_num_dim))
           this%shear_type(:) = 0
    
           IF (ASSOCIATED(this%shear_v0)) THEN
              DEALLOCATE(this%shear_v0)
           END IF
           ALLOCATE(this%shear_v0(d_num_dim,2*d_num_dim))
           this%shear_v0(:,:) = 0.0_MK
           
           IF (ASSOCIATED(this%shear_v)) THEN
              DEALLOCATE(this%shear_v)
           END IF
           ALLOCATE(this%shear_v(d_num_dim,2*d_num_dim))
           this%shear_v(:,:) = 0.0_MK
           
           IF (ASSOCIATED(this%shear_freq)) THEN
              DEALLOCATE(this%shear_freq)
           END IF
           ALLOCATE(this%shear_freq(2*d_num_dim))
           this%shear_freq(:) = 0.0_MK
                      
           IF (ASSOCIATED(this%drag)) THEN
              DEALLOCATE(this%drag)
           END IF
           ALLOCATE(this%drag(d_num_dim,2*d_num_dim))
           this%drag(:,:) = 0.0_MK
           
           this%min_phys(:) = 0.0_MK
           this%max_phys(:) = 0.0_MK
           this%min_phys_t(:) = 0.0_MK
           this%max_phys_t(:) = 0.0_MK

           
        END IF
        
9999    CONTINUE
        
        RETURN       
        
      END SUBROUTINE boundary_set_num_dim
      
      
      SUBROUTINE boundary_set_bcdef(this,d_bcdef,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER, DIMENSION(:)           :: d_bcdef
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        INTEGER                         :: i
        
        stat_info = 0
        
        dim = SIZE(d_bcdef,1)
        
        IF( dim /= 2*this%num_dim) THEN
           PRINT *, "boundary_set_bcdef : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%bcdef(1:dim) = d_bcdef(1:dim)
        
        this%num_peri       = 0
        this%num_sym        = 0
        this%num_wall_sym   = 0
        this%num_wall_solid = 0
        this%num_le         = 0
        
        DO i = 1, dim
           
           SELECT CASE ( d_bcdef(i) )
              
           CASE ( ppm_param_bcdef_periodic ) 
              
              this%num_peri = this%num_peri + 1
              
           CASE ( ppm_param_bcdef_symmetry )
              
              this%num_sym = this%num_sym + 1
              
           CASE( ppm_param_bcdef_wall_sym )
              
              this%num_wall_sym = this%num_wall_sym + 1
              
           CASE( ppm_param_bcdef_wall_solid )
              
              this%num_wall_solid = this%num_wall_solid + 1
              
           CASE( ppm_param_bcdef_LE )
              
              this%num_le = this%num_le + 1
              
           END SELECT ! d_bcdef(i)
           
        END DO ! i = 1,  dim
       
        this%num_shear = this%num_wall_sym  + &
             this%num_wall_solid + &
             this%num_le
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_bcdef
      

      SUBROUTINE boundary_set_shear_type(this,d_shear_type,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER, DIMENSION(:)           :: d_shear_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: i,dim
        
        stat_info = 0
        
        dim = SIZE(d_shear_type,1)
        
        IF( dim /= 2*this%num_dim) THEN
           PRINT *, "boundary_set_shear_type : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_type(1:dim) = d_shear_type(1:dim)
        
        this%num_osci = 0
        
        DO i =1, dim
           IF ( d_shear_type(i) == 2 ) THEN
              this%num_osci = this%num_osci + 1
           END IF
        END DO
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_type
      
      
      SUBROUTINE boundary_set_shear_rate(this,d_shear_rate,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:,:)        :: d_shear_rate
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim1
        INTEGER                         :: dim2

        stat_info = 0
        
        dim1 = SIZE(d_shear_rate,1)
        dim2 = SIZE(d_shear_rate,2)
        
        IF( dim1 /= this%num_dim .OR. &
             dim2 /= this%num_dim ) THEN
           PRINT *, "boundary_set_shear_rate : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_rate(1:dim1,1:dim2) = &
             d_shear_rate(1:dim1,1:dim2) 
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_rate


      SUBROUTINE boundary_set_shear_length(this,d_shear_length,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:,:)        :: d_shear_length
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim1
        INTEGER                         :: dim2

        stat_info = 0
        
        dim1 = SIZE(d_shear_length,1)
        dim2 = SIZE(d_shear_length,2)
        
        IF( dim1 /= this%num_dim .OR. &
             dim2 /= 2*this%num_dim ) THEN
           PRINT *, "boundary_set_shear_length : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_length(1:dim1,1:dim2) = &
             d_shear_length(1:dim1,1:dim2) 
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_length
      
      
      SUBROUTINE boundary_set_shear_v0(this,d_shear_v0,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:,:)        :: d_shear_v0
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim1
        INTEGER                         :: dim2

        stat_info = 0
        
        dim1 = SIZE(d_shear_v0,1)
        dim2 = SIZE(d_shear_v0,2)
        
        IF( dim1 /= this%num_dim .OR. &
             dim2 /= 2*this%num_dim ) THEN
           PRINT *, "boundary_set_shear_v0 : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_v0(1:dim1,1:dim2) = &
             d_shear_v0(1:dim1,1:dim2)
        
        CALL boundary_set_shear_v(this,d_shear_v0(1:dim1,1:dim2),stat_info)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_v0
      
      
      SUBROUTINE boundary_set_shear_v(this,d_shear_v,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:,:)        :: d_shear_v
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim1
        INTEGER                         :: dim2

        stat_info = 0
        
        dim1 = SIZE(d_shear_v,1)
        dim2 = SIZE(d_shear_v,2)
        
        IF( dim1 /= this%num_dim .OR. &
             dim2 /= 2*this%num_dim ) THEN
           PRINT *, "boundary_set_shear_v : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_v(1:dim1,1:dim2) = &
             d_shear_v(1:dim1,1:dim2) 
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_v
      
     
      SUBROUTINE boundary_set_shear_freq(this,d_shear_freq,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:)          :: d_shear_freq
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim

        stat_info = 0
        
        dim = SIZE(d_shear_freq,1)
        
        IF( dim /= 2*this%num_dim ) THEN
           PRINT *, "boundary_set_shear_freq : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%shear_freq(1:dim) = &
             d_shear_freq(1:dim) 
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_shear_freq

      
      SUBROUTINE boundary_set_wall_rho_type(this,d_rho_type,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER                         :: d_rho_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        IF( d_rho_type < 0  .OR. d_rho_type > 1) THEN
           PRINT *, "boundary_set_wall_rho_type : ", &
                "Wrong rho type !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%rho_type =d_rho_type
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_wall_rho_type
      

      SUBROUTINE boundary_set_wall_noslip_type(this,d_noslip_type,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER                         :: d_noslip_type
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%noslip_type =d_noslip_type
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_wall_noslip_type
      
      
      SUBROUTINE boundary_set_dout(this,d_dout,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), INTENT(IN)            :: d_dout
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%dout = d_dout
        
        RETURN
        
      END SUBROUTINE boundary_set_dout
      
      
      SUBROUTINE boundary_set_min_phys(this,d_min_phys,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:)          :: d_min_phys
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_min_phys)
        
        IF( dim /= this%num_dim ) THEN
           PRINT *, "boundary_set_min_phys : ",&
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%min_phys(1:dim) =d_min_phys(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_min_phys

      
      SUBROUTINE boundary_set_max_phys(this,d_max_phys,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:)          :: d_max_phys
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_max_phys)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "boundary_set_max_phys : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%max_phys(1:dim) =d_max_phys(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_max_phys
      

      SUBROUTINE boundary_set_min_phys_t(this,d_min_phys_t,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:)          :: d_min_phys_t
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_min_phys_t)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "boundary_set_min_phys_t : ",&
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%min_phys_t(1:dim) =d_min_phys_t(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_min_phys_t

      
      SUBROUTINE boundary_set_max_phys_t(this,d_max_phys_t,stat_info)
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK), DIMENSION(:)          :: d_max_phys_t
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        dim = SIZE(d_max_phys_t)
        
        IF( dim /= this%num_dim) THEN
           PRINT *, "boundary_set_max_phys_t : ", &
                "Wrong dimension !"
           stat_info = -1
           GOTO 9999
        END IF
        
        this%max_phys_t(1:dim) =d_max_phys_t(1:dim)
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE  boundary_set_max_phys_t
      
      
      SUBROUTINE boundary_set_num_part_wall_solid(this,num_part,stat_info)
        
        !---------------------------------------
        ! Set the num of wall boundary particles
        ! created by MCF.
        !---------------------------------------
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        INTEGER, INTENT(IN)             :: num_part
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        this%num_part_wall_solid = num_part
        
        RETURN
        
      END SUBROUTINE  boundary_set_num_part_wall_solid
      
