      SUBROUTINE io_write_boundary(this,&
           rank,step,time,d_boundary,stat_info)
        !----------------------------------------------------
        !  Subroutine   :  io_write_boundary
        !----------------------------------------------------
        !
        !  Purpose      :  Write information into output files.
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : 
        !                 V0.1 13.01.2009, original version,
        !
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
	! Arguments
	!----------------------------------------------------
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)             :: rank
        INTEGER,  INTENT(IN)	        :: step
        REAL(MK), INTENT(IN)	        :: time
        TYPE(Boundary), INTENT(IN)      :: d_boundary
        INTEGER,  INTENT(OUT)	        :: stat_info
        
        !----------------------------------------------------
        !  Local Variables
        !----------------------------------------------------
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: Brownian
        INTEGER                         :: num_dim
        INTEGER                         :: num_shear
        INTEGER                         :: num_wall_solid
        INTEGER                         :: num_wall_sym
        INTEGER,DIMENSION(:), POINTER   :: bcdef
        REAL(MK),DIMENSION(:,:),POINTER :: drag

        REAL(MK),DIMENSION(:,:),POINTER :: drag_p
        REAL(MK),DIMENSION(:,:),POINTER :: drag_v
        REAL(MK),DIMENSION(:,:),POINTER :: drag_r
        
        REAL(MK),DIMENSION(:), POINTER  :: output
        INTEGER                         :: num
        INTEGER                         :: j_start,j_end
        CHARACTER(len=MAX_CHAR)	        :: cbuf,fbuf
        INTEGER			        :: i

        !----------------------------------------------------
        ! Initialization.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        Brownian = &
             control_get_Brownian(this%ctrl,stat_info_sub)
        
        NULLIFY(bcdef)
        NULLIFY(drag)        
        
        NULLIFY(drag_p)
        NULLIFY(drag_v)
        NULLIFY(drag_r)

        NULLIFY(output)
        
        IF ( rank /= 0 ) THEN
           PRINT *, &
                "io_write_boundary can only be used by root process !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Get 
        ! number of dimension;
        ! number of shear;
        ! boundary conditions;
        ! drags on the walls.
        !----------------------------------------------------
        
        num_dim   = boundary_get_num_dim(d_boundary,stat_info_sub)
        num_shear = boundary_get_num_shear(d_boundary,stat_info_sub)

        !----------------------------------------------------
        ! If there is no wall, nothing to write yet.
        !----------------------------------------------------
        
        IF ( num_shear > 0) THEN
           
           num_wall_solid = &
                boundary_get_num_wall_solid(d_boundary,stat_info_sub)
           num_wall_sym = &
                boundary_get_num_wall_sym(d_boundary,stat_info_sub)
           
           num = 1 + &
                num_dim*num_wall_solid + &
                num_dim*num_wall_sym
           
           CALL boundary_get_bcdef(d_boundary,bcdef,stat_info_sub)
           CALL boundary_get_drag(d_boundary,drag,stat_info_sub)
           
#ifdef __IO_WALL_FORCE_SEPARATE
           CALL boundary_get_drag_p(d_boundary,drag_p,stat_info_sub)
           CALL boundary_get_drag_v(d_boundary,drag_v,stat_info_sub)
           num = num + &
                2*num_dim*num_wall_solid + &
                2*num_dim*num_wall_sym
           IF ( Brownian ) THEN
              CALL boundary_get_drag_r(d_boundary,drag_r,stat_info_sub)
              num = num + &
                   num_dim*num_wall_solid + &
                   num_dim*num_wall_sym
           END IF
#endif 
           !-------------------------------------------------
           ! Allocate memory for output data.
           !-------------------------------------------------
           
           ALLOCATE(output(num))
           
           !-------------------------------------------------
           ! Record current time.
           !-------------------------------------------------
           
           j_start = 1
           j_end   = 1
           
           output(j_start) = time
           
           DO i = 1, 2*num_dim
              
              IF ( bcdef(i) == ppm_param_bcdef_wall_sym .OR. &
                   bcdef(i) == ppm_param_bcdef_wall_solid ) THEN 
                 
                 !-------------------------------------------
                 ! Record force/drag.
                 !-------------------------------------------
                 
                 j_start = j_end + 1
                 j_end   = j_start + num_dim -1
                 
                 output(j_start:j_end) = &
                      drag(1:num_dim,i)
                 
              END IF ! bcdef
              
           END DO
           
#ifdef __IO_WALL_FORCE_SEPARATE
           
           DO i = 1, 2*num_dim
              
              IF ( bcdef(i) == ppm_param_bcdef_wall_sym .OR. &
                   bcdef(i) == ppm_param_bcdef_wall_solid ) THEN 
                 
                 !-------------------------------------------
                 ! Record force/drag separately.
                 !-------------------------------------------
                 
                 j_start = j_end + 1
                 j_end   = j_start + num_dim -1
                 
                 output(j_start:j_end) = &
                      drag_p(1:num_dim,i)
                 
                 j_start = j_end + 1
                 j_end   = j_start + num_dim -1
                 
                 output(j_start:j_end) = &
                      drag_v(1:num_dim,i)

                 IF ( Brownian ) THEN
                    
                    j_start = j_end + 1
                    j_end   = j_start + num_dim -1
                    
                    output(j_start:j_end) = &
                         drag_r(1:num_dim,i)
                    
                 END IF

                 
              END IF ! bcdef
              
           END DO
           
#endif
           
           WRITE(fbuf, '(A,I2,A)' ) , '(I9,', j_end, 'E16.8)'
           WRITE(cbuf,fbuf )  step, output(1:j_end)
           
           WRITE(UNIT=this%boundary_unit,FMT='(A)',&
                IOSTAT=stat_info_sub)  TRIM(cbuf)
           
           IF( stat_info_sub /= 0 ) THEN
              PRINT *,"io_write_boundary : ",&
                   "Writting into boundary file failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF

        IF(ASSOCIATED(output)) THEN
           DEALLOCATE(output)
        END IF
        
        IF(ASSOCIATED(drag)) THEN
           DEALLOCATE(drag)
        END IF
        
#ifdef __IO_WALL_FORCE_SEPARATE
        IF(ASSOCIATED(drag_p)) THEN
           DEALLOCATE(drag_p)
        END IF
        IF(ASSOCIATED(drag_v)) THEN
           DEALLOCATE(drag_v)
        END IF
        IF(ASSOCIATED(drag_r)) THEN
           DEALLOCATE(drag_r)
        END IF
#endif
        
        RETURN 
        
      END SUBROUTINE io_write_boundary
      
