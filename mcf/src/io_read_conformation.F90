      SUBROUTINE io_read_conformation(this,d_rank,d_particles,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_read_conformation
        ! 
        ! Purpose     : Reading particles' conformation
        !               tensor from file.
        !
        ! Revision    : V0.2 04.12 2009, check the work flow
        !               and supply with more comments.
        !
        !               V0.1 05.08 2009, original version.
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
        
        TYPE(IO), INTENT(IN)                    :: this
        INTEGER, INTENT(IN)                     :: d_rank
        TYPE(Particles), INTENT(INOUT)          :: d_particles
        INTEGER, INTENT(OUT)                    :: stat_info
        
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        REAL(MK), DIMENSION(:,:), POINTER       :: ct
        INTEGER                                 :: ilenread
        CHARACTER(MAX_CHAR)                     :: cbuf
        LOGICAL                                 :: lExist
        TYPE(Physics),POINTER                   :: phys
        INTEGER                                 :: num_dim
        INTEGER                                 :: num_part
        INTEGER                                 :: iline

#ifdef __DEBUG
        !----------------------
        !  Debug variables.
        !----------------------
        INTEGER                                 :: debug_flag
        INTEGER                                 :: debug_threshold
        REAL(MK)                                :: time_routine_start
#endif       

        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(ct)
        
        
        IF( d_rank /= 0) THEN
           PRINT *, "io_read_conformation : " , &
                "can only be called by root process ! "
           stat_info  = -1
           GOTO 9999
        END IF
        
        
#ifdef __DEBUG
        !----------------------
        !  Debug purpose.
        !----------------------
        debug_threshold = 1        
        debug_flag = debug_get_flag(global_debug,stat_info_sub)
        IF(debug_flag > 1 .OR. debug_flag > debug_threshold)  THEN
           CALL debug_substart(global_debug,&
                d_rank,'io_read_conformation',&
                time_routine_start,stat_info_sub)
        END IF
#endif
        
        CALL particles_get_phys(d_particles,phys,stat_info_sub)        
        num_dim = physics_get_num_dim(phys,stat_info_sub)
        
	!----------------------------------------------------
      	! Check if name of particle conformation file is empty.
      	!----------------------------------------------------
        
        ilenread = LEN_TRIM(this%read_conformation_file)
        
        IF ( ilenread < 1 ) THEN
           PRINT *,'io_read_conformation : ',&
                'No file name given !'
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
      	! Check if the particle conformation file exists. 
      	!----------------------------------------------------
        
        INQUIRE(FILE=this%read_conformation_file,EXIST=lExist)
	
        IF (.NOT.lExist) THEN
           WRITE(cbuf,'(2A)')'No such file: ', &
                this%read_conformation_file(1:ilenread)
           PRINT *, 'io_read_conformation : "', cbuf
           stat_info = -1
           GOTO 9999
        END IF
        
      	!----------------------------------------------------
      	! Open the file.
      	!----------------------------------------------------
        
        OPEN(this%read_conformation_unit,&
             FILE=this%read_conformation_file,&
             IOSTAT=stat_info_sub,ACTION='READ')
        
        IF (stat_info_sub /= 0) THEN
           WRITE(cbuf,'(2A)')'Failed to open file: ',&
                this%read_conformation_file(1:ilenread)
           PRINT *,'io_read_conformation : ', cbuf
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Scan the first file line for conformaiotn number.
        !----------------------------------------------------
        
        READ(this%read_conformation_unit,*,END=9999,ERR=200) &
             num_part
        
        !----------------------------------------------------
        ! Allocate memory according to conformation number.
        !----------------------------------------------------
        
        IF( num_part > 0 ) THEN
           
           ALLOCATE(ct(num_dim**2,num_part))
           
        END IF
        
        !----------------------------------------------------
        ! Scan the file line by line and read conformation.
        !----------------------------------------------------
        
        iline = 0
        
        DO
           
           !-------------------------------------------------
           ! Increase line counter.
           !-------------------------------------------------
           
           iline = iline + 1
           
           
           IF(iline > num_part ) THEN
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Read information of current line.
           !-------------------------------------------------
           
           READ(this%read_conformation_unit,*,END=9999,ERR=200) &
                ct(1:num_dim**2,iline)
           
        END DO
        
200     CONTINUE
        
        PRINT *, "io_read_conformation : ", &
             "reading conformation tensor has problem at line ",&
             iline
        stat_info = -1
        GOTO 9999
        
9999    CONTINUE
        
        IF (iline /= num_part+1) THEN
           
           PRINT *, "io_read_conformation : ", &
                "actual number of lines is not equal to number given !"
           stat_info = -1
           
        ELSE
           
           CALL particles_set_ct(d_particles,ct,num_part,stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Release dynamic memory.
        !----------------------------------------------------
        
        IF(ASSOCIATED(ct)) THEN
           DEALLOCATE(ct)
        END IF
        
        
#ifdef __DEBUG        
        IF(debug_flag > 1 .OR. debug_flag > debug_threshold)  THEN
           CALL debug_substop(global_debug,d_rank,&
                'io_read_conformation',&
                time_routine_start,stat_info_sub)
        END IF
#endif        
        
        RETURN
        
      END SUBROUTINE io_read_conformation
