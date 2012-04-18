      SUBROUTINE io_write_colloid_separate(this,&
           rank,step,time,d_colloid,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_write_colloid_separate
        !----------------------------------------------------
        !
        ! Purpose     : Write information into output files.
        !
        ! Reference   :
        !
        ! Remark      : 
        !
        ! Revisions   : V0.1 13.01.2009, original version,
        !               each colloid each file, done by
        !               rank = 0
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
        
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)             :: rank
        INTEGER,  INTENT(IN)	        :: step
        REAL(MK), INTENT(IN)	        :: time
        TYPE(Colloid), INTENT(IN)       :: d_colloid
        INTEGER,  INTENT(OUT)	        :: stat_info
        
        !----------------------------------------------------
        !  Local Variables
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: read_external
        INTEGER                         :: num_dim
        INTEGER                         :: num_colloid
        LOGICAL                         :: translate
        LOGICAL                         :: rotate
        REAL(MK),DIMENSION(:,:),POINTER :: x
        REAL(MK),DIMENSION(:,:,:),POINTER:: v
#ifdef __DRAG_PART
        REAL(MK),DIMENSION(:,:),POINTER :: drag_lub
        REAL(MK),DIMENSION(:,:),POINTER :: drag_repul
#endif
        REAL(MK),DIMENSION(:,:),POINTER :: drag
        REAL(MK),DIMENSION(:,:,:),POINTER::rot_matrix
        REAL(MK),DIMENSION(:,:),POINTER :: theta
        REAL(MK),DIMENSION(:,:,:),POINTER:: omega
        REAL(MK),DIMENSION(:,:),POINTER :: torque
        
        INTEGER                         :: data_dim
        REAL(MK),DIMENSION(:,:),POINTER :: output
        INTEGER                         :: j_start,j_end
        CHARACTER(len=2*MAX_CHAR)       :: cbuf,fbuf
        INTEGER			        :: i, j
        
        !---------------------------------------------------
        ! Initialization.
        !----------------------------------------------------
        
	stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(x)
        NULLIFY(v)
        NULLIFY(rot_matrix)
        NULLIFY(theta)
        NULLIFY(omega)
#ifdef __DRAG_PART
        NULLIFY(drag_lub)
        NULLIFY(drag_repul)        
#endif
        NULLIFY(drag)
        NULLIFY(torque)
        NULLIFY(output)


	!----------------------------------------------------
        ! Can not be written by non 0 rank.
	!----------------------------------------------------
        
        IF ( rank /= 0 ) THEN
           PRINT *, &
                "io_write_colloid can only be used by root process!"
           stat_info = -1
           GOTO 9999           
        END IF
        
	!----------------------------------------------------
        ! Get control and colloidal quantity parameters.
	!----------------------------------------------------
        
        read_external   = &
             control_get_read_external(this%ctrl,stat_info_sub)
        
        num_dim     = &
             colloid_get_num_dim(d_colloid,stat_info_sub)
        num_colloid = &
             colloid_get_num_colloid(d_colloid,stat_info_sub)
        translate = &
             colloid_get_translate(d_colloid,stat_info_sub)
        rotate    = &
             colloid_get_rotate(d_colloid,stat_info_sub)
        CALL colloid_get_x(d_colloid,x,stat_info_sub)
        CALL colloid_get_v(d_colloid,v,stat_info_sub)
#ifdef __DRAG_PART
        CALL colloid_get_drag_lub(d_colloid,drag_lub,stat_info_sub)
        CALL colloid_get_drag_repul(d_colloid,drag_repul,stat_info_sub)
#endif
        CALL colloid_get_drag(d_colloid,drag,stat_info_sub)

        IF ( num_dim == 3) THEN
           CALL colloid_get_accumulation_matrix(d_colloid,&
                rot_matrix,stat_info_sub)
        END IF
        CALL colloid_get_theta(d_colloid,theta,stat_info_sub)
        CALL colloid_get_omega(d_colloid,omega,stat_info_sub)
        CALL colloid_get_torque(d_colloid,torque,stat_info_sub)
        
        
        !----------------------------------------------------
        ! Allocate memory for output collective data.
        ! x, v, F, theta, omega, torque
        !----------------------------------------------------
        
        IF ( num_dim == 2 ) THEN
           data_dim = 15
        ELSE IF (num_dim ==3 ) THEN
           data_dim = 32
        END IF
        
        ALLOCATE(output(data_dim,num_colloid),STAT=stat_info_sub)
        
        IF (stat_info_sub /=0) THEN          
           PRINT *, 'io_write_colloid_separate : ',&
                'Allocating buffer for output failed !'
           stat_info = -1
           GOTO 9999          
        END IF
        
        !----------------------------------------------------
        ! Set starting indices for recording data
        !----------------------------------------------------
        
        j_start = 0
        j_end   = 0
        
        output(1,1:num_colloid) = time
        
        j_start = 1
        j_end   = 1
        !----------------------------------------------------
        ! Record position and velocity for free translating
        ! colloids, i.e., it they do not move, not need
        ! to record x, or v at all.
        !----------------------------------------------------
        
        IF ( translate ) THEN
           
           j_start = j_end + 1
           j_end   = j_start + num_dim -1
           
           output(j_start:j_end,1:num_colloid) = &
                x(1:num_dim, 1:num_colloid)
           
           j_start = j_end + 1
           j_end   = j_start + num_dim -1
              
           output(j_start:j_end,1:num_colloid) = &
                v(1:num_dim,1:num_colloid,1)
           
        END IF ! translate
        
        
        !----------------------------------------------------
        ! Record drag no matter whether it moves.
        !----------------------------------------------------
        
#if __DRAG_PART
        !----------------------------------------------------
        ! If different contribution to the total drag
        ! is interesting, then record them also.
        ! 1: Lubrication correction
        ! 2: repulsive force.
        !----------------------------------------------------
        
        j_start = j_end + 1
        j_end   = j_start + num_dim -1
        
        output(j_start:j_end,1:num_colloid) = &
             drag_lub(1:num_dim,1:num_colloid)
        
        j_start = j_end + 1
        j_end   = j_start + num_dim -1
        
        output(j_start:j_end,1:num_colloid) = &
             drag_repul(1:num_dim,1:num_colloid)
#endif
        
        j_start = j_end + 1
        j_end   = j_start + num_dim -1
        
        output(j_start:j_end,1:num_colloid) = &
             drag(1:num_dim,1:num_colloid)
        !----------------------------------------------------
        ! If colloid rotates, record rotated angle(2D) 
        ! or rotation matrix(3D), and angular velocity.
        ! No need if they do not rotate.
        !----------------------------------------------------
           
        IF ( rotate ) THEN 
              
           IF ( num_dim == 2 ) THEN
                 
              !----------------------------------------------
              ! For 2D, 
              ! rotated angle is 1D.
              ! angular velocity is 1D also.
              !----------------------------------------------
              
              j_start = j_end + 1
              j_end   = j_start
              
              output(j_end,1:num_colloid) = &
                   theta(3,1:num_colloid)
              
              j_start = j_end + 1
              j_end   = j_start
              
              output(j_end,1:num_colloid) = &
                   omega(3,1:num_colloid,1)
              
           ELSE IF ( num_dim == 3 ) THEN
              
              !----------------------------------------------
              ! For 3D, 
              ! record theta
              ! record angular velocity.
              !----------------------------------------------
              
              j_start = j_end + 1
              j_end   =  j_start + 8
              
              DO i = 1, 3
                 DO j =1, 3
                    output(j_start+(i-1)*3+j-1,1:num_colloid) = &
                         rot_matrix(j,i,1:num_colloid)
                 END DO
              END DO
              
              !j_start = j_end + 1
              !j_end   = j_start + num_dim -1
              
              !output(j_start:j_end,1:num_colloid) = &
              !     theta(1:3,1:num_colloid)
              
              j_start = j_end + 1
              j_end   = j_start + num_dim - 1
              
              output(j_start:j_end,1:num_colloid) = &
                   omega(1:3,1:num_colloid,1)
              
           END IF ! num_dim
           
        END IF ! rotate
        
        
        !----------------------------------------------------
        ! Record torque anyway
        !----------------------------------------------------
        
        IF ( num_dim == 2 ) THEN
           
           !-------------------------------------------------
           ! For 2D, record 1D torque.
           !--------------------------------------------------
           
           j_start = j_end + 1
           j_end   = j_start
           
           output(j_end,1:num_colloid) = &
                torque(3,1:num_colloid)
           
        ELSE
           
           !-------------------------------------------------
           ! For 3D, record 3D torque.
           !-------------------------------------------------
           
           j_start = j_end + 1
           j_end   = j_start + num_dim -1
           
           output(j_start:j_end,1:num_colloid) = &
                torque(1:num_dim,1:num_colloid)
           
        END IF ! num_dim
        

        !----------------------------------------------------
        ! Open file and write data into file.
        !----------------------------------------------------

        WRITE(fbuf, '(A4,I2,A6)' ) , '(I9,', j_end, 'E16.8)'


        DO i = 1, num_colloid
           

           WRITE(cbuf,fbuf) step, output(1:j_end,i)
           
           WRITE(UNIT=this%colloid_unit+i,FMT='(A)',&
                IOSTAT=stat_info_sub)  TRIM(cbuf)
           
           IF( stat_info_sub /= 0 ) THEN
              PRINT *,"io_write_colloid : ",&
                   "Writting into colloid file failed!"
              stat_info = -1
              GOTO 9999
           END IF
           
        END DO ! i = 1, num_colloid


9999    CONTINUE
        
        IF(ASSOCIATED(x)) THEN
           DEALLOCATE(x)
        END IF
        
        IF(ASSOCIATED(v)) THEN
           DEALLOCATE(v)
        END IF
        
        IF(ASSOCIATED(rot_matrix)) THEN
           DEALLOCATE(rot_matrix)
        END IF
        
        IF(ASSOCIATED(theta)) THEN
           DEALLOCATE(theta)
        END IF

        IF(ASSOCIATED(omega)) THEN
           DEALLOCATE(omega)
        END IF
        
#if __DRAG_PART
        IF(ASSOCIATED(drag_lub)) THEN
           DEALLOCATE(drag_lub)
        END IF
        IF(ASSOCIATED(drag_repul)) THEN
           DEALLOCATE(drag_repul)
        END IF
#endif

        IF(ASSOCIATED(drag)) THEN
           DEALLOCATE(drag)
        END IF

        IF(ASSOCIATED(torque)) THEN
           DEALLOCATE(torque)
        END IF
        
        IF(ASSOCIATED(output)) THEN
           DEALLOCATE(output)
        END IF
        
        RETURN 
        
      END SUBROUTINE io_write_colloid_separate
      
