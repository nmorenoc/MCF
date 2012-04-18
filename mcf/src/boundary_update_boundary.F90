      SUBROUTINE boundary_update_boundary(this,time_r,time_v,stat_info)
        !----------------------------------------------------
        ! Subroutine  :  boundary_update_boundary
        !----------------------------------------------------
        !
        ! Purpose     :  For solid/symmetry wall, update 
        !                its velocity according to its type;
        !                For Lees-Edwards boundary, update
        !                its velocity and sheared length 
        !                on each side.
        !
        ! Reference   :
        !
        ! Remark      : For the moment, we use always imposed
        !               shear rate set-up, instead of imposed
        !               shear stress on the wall, therefore,
        !               velocity of the wall is known already 
        !               according to time and there is no need 
        !               to integrate wall position or velocity.
        !
        ! Revisions   : V0.2 18.11.2009, including
        !               Lees-Edwards boundary.
        !                 
        !               V0.1 13.10.2009, original version.
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
        
        TYPE(Boundary), INTENT(INOUT)   :: this
        REAL(MK),INTENT(IN)             :: time_r
        REAL(MK),INTENT(IN)             :: time_v
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------

        INTEGER                         :: i, dim
        INTEGER                         :: num_wall
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization.
        !----------------------------------------------------

        stat_info     = 0
        stat_info_sub = 0
        
        dim       = this%num_dim
        num_wall  = &
             boundary_get_num_wall(this,stat_info_sub)
        
        !----------------------------------------------------
        ! Check current time given reasonable.
        !----------------------------------------------------
        
        IF( time_r < 0.0_MK .OR. time_v < 0.0_MK ) THEN
           PRINT *, "boundary_update_boundary: ", &
                "Negative time is impossible !"
           !stat_info = -1
           !GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! If there is shear boundary
        !----------------------------------------------------
        
        IF ( this%num_shear > 0 ) THEN
        
           !-------------------------------------------------
           ! For oscillating  wall, update shear velocity.
           !-------------------------------------------------
           
           IF ( num_wall > 0 .AND. &
                this%num_osci > 0 )THEN
              
              DO i = 1, 2*dim
                 
                 IF( this%shear_type(i) == 2 ) THEN
                    
                    this%shear_v(1:dim,i) = &
                         this%shear_v0(1:dim,i) * &
                         COS(this%shear_freq(i) * time_v)
                    
                 END IF ! shear_type
                 
              END DO ! i
           
           END IF ! oscillating
        
           !-------------------------------------------------
           ! For Lees-Edwards boundary, update shear length.
           !-------------------------------------------------
           
           IF ( this%num_le > 0 ) THEN
              
              this%shear_length(1:dim,1:2*dim) = 0.0_MK
              
              DO i =1, dim
                 
                 IF( this%bcdef(2*i-1) == ppm_param_bcdef_LE .AND. &
                      this%bcdef(2*i) == ppm_param_bcdef_LE ) THEN
                    
                    this%shear_length(1:dim,2*i-1) = &
                         ( this%shear_v0(1:dim,2*i-1)-&
                         this%shear_v0(1:dim,2*i))* time_r
                    
                    this%shear_length(1:dim,2*i) = &
                         ( this%shear_v0(1:dim,2*i)-&
                         this%shear_v0(1:dim,2*i-1))* time_r
                    
                 END IF ! bcdef
                 
              END DO ! i
              
           END IF ! num_le > 0

        END IF ! num_shear
        
        
9999    CONTINUE
        
        
        RETURN
        
      END SUBROUTINE  boundary_update_boundary
      
