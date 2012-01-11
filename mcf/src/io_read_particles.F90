      SUBROUTINE io_read_particles(this,d_rank,d_particles,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_read_particles
        !----------------------------------------------------
        !
        ! Purpose     : Reading particles' configuration,
        !               i.e. position, velocity, 
        !               mass/number density, mass,
        !               particle ID, species ID from file.
        !
        ! Revision    : V0.2 04.12 2009, check the workflow
        !               and supply with more comments.
        !
        !               V0.1 01.04 2009, original version.
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
        REAL(MK), DIMENSION(:,:), POINTER       :: x
        REAL(MK), DIMENSION(:,:), POINTER       :: v
        REAL(MK), DIMENSION(:), POINTER         :: rho
        REAL(MK), DIMENSION(:), POINTER         :: m
        INTEGER, DIMENSION(:,:), POINTER        :: id
        REAL(MK), DIMENSION(:), POINTER         :: inter_id
        
        INTEGER                                 :: ilenread
        CHARACTER(MAX_CHAR)                     :: cbuf
        LOGICAL                                 :: lExist
        TYPE(Physics),POINTER                   :: phys
        INTEGER                                 :: num_dim
        INTEGER                                 :: num_id
        INTEGER                                 :: num_part
        INTEGER                                 :: iline
        
#ifdef __DEBUG
        !----------------------------------------------------
        ! Debug variables.
        !----------------------------------------------------

        INTEGER                                 :: debug_flag
        INTEGER                                 :: debug_threshold
        REAL(MK)                                :: time_routine_start
#endif       
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        
        NULLIFY(x)
        NULLIFY(v)
        NULLIFY(rho)
        NULLIFY(m)
        NULLIFY(id)
        NULLIFY(inter_id)
        NULLIFY(phys)
        
        
        IF( d_rank /= 0 ) THEN
           PRINT *, "io_read_particles : " , &
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
                d_rank,'io_read_particles',&
                time_routine_start,stat_info_sub)
        END IF
#endif
        
        
	!----------------------------------------------------
        ! num_dim : number of dimension
        ! num_id  : number of different type of IDs for particle.
	!----------------------------------------------------

        num_dim = particles_get_num_dim(d_particles,stat_info_sub)
        num_id  = particles_get_num_id(d_particles,stat_info_sub)
        
	!----------------------------------------------------
      	! Check if name of particle file is empty.
      	!----------------------------------------------------
        
        ilenread = LEN_TRIM(this%read_particles_file)
        
        IF ( ilenread < 1 ) THEN
           PRINT *,'io_read_particles : ',&
                'No file name is given !'
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
      	! Check if the particle file exists. 
      	!----------------------------------------------------
        
        INQUIRE(FILE=this%read_particles_file,EXIST=lExist)
	
        IF (.NOT.lExist) THEN
           WRITE(cbuf,'(2A)')'No such file: ', &
                this%read_particles_file(1:ilenread)
           PRINT *, 'io_read_particles : "', cbuf
           stat_info = -1
           GOTO 9999
        END IF
        
      	!----------------------------------------------------
      	! Open the file.
      	!----------------------------------------------------
        
        OPEN(this%read_particles_unit,&
             FILE=this%read_particles_file,&
             IOSTAT=stat_info_sub,ACTION='READ')
        
        IF (stat_info_sub /= 0) THEN
           WRITE(cbuf,'(2A)')'Failed to open file: ',&
                this%read_particles_file(1:ilenread)
           PRINT *,'io_read_particles : ', cbuf
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Scan the first file line for particle number.
        !----------------------------------------------------
        
        READ(this%read_particles_unit,*,END=9999,ERR=200) &
             num_part
        
        !----------------------------------------------------
        ! Allocate memory according to particle number.
        !----------------------------------------------------
        
        IF( num_part > 0) THEN
           
           ALLOCATE(x(num_dim,num_part))
           ALLOCATE(v(num_dim,num_part))
           ALLOCATE(rho(num_part))
           ALLOCATE(m(num_part))
           ALLOCATE(id(num_id,num_part))
           ALLOCATE(inter_id(num_id))
           
        END IF
        
        !----------------------------------------------------
        ! Scan the file line by line and read particles.
        !----------------------------------------------------
        
        iline = 0
        
        DO
           !-------------------------------------------------
           ! Increase line counter.
           !-------------------------------------------------
           
           iline = iline + 1
           
           
           IF( iline > num_part ) THEN
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Read information of current line.
           !-------------------------------------------------
           
           READ(this%read_particles_unit,*,END=9999,ERR=200) &
                x(1:num_dim,iline), v(1:num_dim,iline),&
                rho(iline), m(iline), inter_id(1:num_id)
           
           !-------------------------------------------------
           ! Convert floating point number of ID to interger.
           !-------------------------------------------------
           
           id(1:num_id,iline) = INT(inter_id(1:num_id))
           
        END DO
        
200     CONTINUE
        
        PRINT *, "io_read_particles : ", &
             "reading particles has problem at line ",&
             iline
        stat_info = -1
        GOTO 9999
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! If actual read lines is not equal to the number
        ! of particles given in the beginning of file,
        ! then there is something some.
        ! Otherwise, save the information of particles.
        !----------------------------------------------------

        IF ( iline /= num_part+1 ) THEN
           
           PRINT *, "io_read_particles : ", &
                "actual number of lines is not equal to number given !"
           stat_info = -1
           
        ELSE
           
           CALL particles_init_global_exter(d_particles,&
                x,v,rho,m,id,num_part,stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Release dynamic memory.
        !----------------------------------------------------
        
        IF(ASSOCIATED(x)) THEN
           DEALLOCATE(x)
        END IF
        
        IF(ASSOCIATED(v)) THEN
           DEALLOCATE(v)
        END IF
        
        IF(ASSOCIATED(rho)) THEN
           DEALLOCATE(rho)
        END IF
        
        IF(ASSOCIATED(m)) THEN
           DEALLOCATE(m)
        END IF
        
        IF(ASSOCIATED(id)) THEN
           DEALLOCATE(id)
        END IF
        
        IF(ASSOCIATED(inter_id)) THEN
           DEALLOCATE(inter_id)
        END IF
        
        
#ifdef __DEBUG        
        IF(debug_flag > 1 .OR. debug_flag > debug_threshold)  THEN
           CALL debug_substop(global_debug,d_rank,&
                'io_read_particles',&
                time_routine_start,stat_info_sub)
        END IF
#endif        
        
        RETURN
        
      END SUBROUTINE io_read_particles
      
