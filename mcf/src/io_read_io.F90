      SUBROUTINE io_read_io_config(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : io_read_io_config
        !----------------------------------------------------
        !
        ! Purpose     :  Reading io configuration,
        !                file names, format, frequncy, etc.
        !
        !  Revision   : V0.1 01.04.2009
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
        ! Argument
        !----------------------------------------------------
        
        TYPE(IO), INTENT(INOUT)            :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables 
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                        	:: i,j,idx
        INTEGER                        	:: ilen,ios,iline,ilenio
        CHARACTER(LEN=MAX_CHAR)         :: cbuf,cvalue,carg
        LOGICAL                        	:: lExist
        
        
#ifdef __DEBUG
        INTEGER                         :: debug_threshold
#endif
        
	!----------------------------------------------------
      	! Initialization.
      	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
#ifdef __DEBUG
        debug_threshold = 0
        IF ( debug_threshold > 1) THEN           
           PRINT *, "io_read_io_config : starting "
        END IF
#endif
        
        !----------------------------------------------------
        ! Check if the name of io config file exists
      	!----------------------------------------------------
        
        ilenio = LEN_TRIM(this%io_config_file)
        IF (ilenio < 1) THEN
           PRINT *,'io_read_io_config : ',&
                'No io config file name given!'
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check if the io config file exists 
        !----------------------------------------------------
        
        INQUIRE(FILE=this%io_config_file,EXIST=lExist)
        IF (.NOT.lExist) THEN
           WRITE(cbuf,'(2A)')'No such file: ',&
                this%io_config_file(1:ilenio)
           PRINT *, 'io_read_io_config : "', cbuf
           stat_info = -1
           GOTO 200
        END IF
        
        !----------------------------------------------------
        ! Open the file
        !----------------------------------------------------
        
        OPEN(this%io_config_unit,FILE=this%io_config_file,&
             IOSTAT=ios,ACTION='READ')
        
        IF (ios /= 0) THEN
           WRITE(cbuf,'(2A)')'Failed to open file: ',&
                this%io_config_file(1:ilenio)
           PRINT *, 'io_read_io_config : ', cbuf
           stat_info = -1
           GOTO 200
        END IF
        
        !----------------------------------------------------
        ! Scan the file line by line
        !----------------------------------------------------
        
        iline = 0
        
        DO
           !-------------------------------------------------
           ! Increase line counter
           !-------------------------------------------------
           
           iline = iline + 1
           
           !-------------------------------------------------
           ! Read current line
           !-------------------------------------------------
           
           READ(this%io_config_unit,'(A)',END=9999,ERR=200) cbuf 
           ilen = LEN_TRIM(cbuf)
           
           !-------------------------------------------------
           ! Extensive Debug: print each line of the file
           !-------------------------------------------------
           
#ifdef __DEBUG
           IF ( debug_threshold > 3 ) THEN
              PRINT *,'io_read_io_config : ', cbuf(1:ilen)
           END IF
#endif          
           
           !-------------------------------------------------
           ! Skip empty line or comment line with symbol
           ! at first place 
           !-------------------------------------------------
           
           IF (ilen < 1 .OR. cbuf(1:1) == '#' ) THEN
              CYCLE
           END IF
           
           !-------------------------------------------------
           ! Remove spaces of the line being read
           !-------------------------------------------------
          
           j = 0
           DO i=1,ilen
              IF (cbuf(i:i) /=' ' .AND. cbuf(i:i) /= '\t' ) THEN
                 j = j + 1
                 cbuf(j:j) = cbuf(i:i)
              END IF
           END DO
           
           !-------------------------------------------------
           ! Update length of string
           !-------------------------------------------------
           
           ilen = j
           
           !-------------------------------------------------
           ! Skip comment line with symbol not at first place 
           !-------------------------------------------------
           
           IF (cbuf(1:1) == '#' ) THEN
              CYCLE
           END IF
           
           !-------------------------------------------------
           ! Find position of '='
           !-------------------------------------------------
           
           idx = INDEX(cbuf,'=')
           
           !-------------------------------------------------
           ! Exit if '=' is missing
           !-------------------------------------------------
           
           IF (idx < 0) THEN
              WRITE(cbuf,'(A,I5)')'Incorrect line: ',iline
              Print *,'io_read_io_config : ', cbuf
              stat_info = -1
              GOTO 200
           END IF
           
           !-------------------------------------------------
           ! Get argument and value
           !-------------------------------------------------
           
           carg   = ADJUSTL(cbuf(1:idx-1))
           cvalue = ADJUSTL(cbuf(idx+1:ilen))
           
           !-------------------------------------------------
           ! Convert to upper case
           !-------------------------------------------------
           
           CALL tool_uppercase(this%io_tool,carg,idx-1,stat_info)
           
#ifdef __DEBUG
           IF ( debug_threshold > 2) THEN
              PRINT *, 'io_read_io_config : ', carg
              PRINT *, 'io_read_io_config : ', cvalue
           END IF
#endif     
           
           !-------------------------------------------------
           ! Now read in the parameters.
           !-------------------------------------------------
           
           IF (carg == 'READ_PARTICLES_FILE') THEN
              
              this%read_particles_file = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'READ_PARTICLES_FMT') THEN
	      
              this%read_particles_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'READ_CONFORMATION_FILE') THEN
              
              this%read_conformation_file = TRIM(ADJUSTL(cvalue))
             
           ELSE IF (carg == 'READ_CONFORMATION_FMT') THEN
	      
              this%read_conformation_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE  IF (carg == 'OUTPUT_PARTICLES_RELAX_FILE') THEN
              
              this%output_particles_relax_file = TRIM(ADJUSTL(cvalue))
             
           ELSE IF (carg == 'OUTPUT_PARTICLES_RELAX_FMT') THEN
              
              this%output_particles_relax_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'OUTPUT_PARTICLES_RELAX_FREQ_STEP') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                   this%output_particles_relax_freq_step
              
           ELSE  IF (carg == 'OUTPUT_PARTICLES_FILE') THEN
              
              this%output_particles_file = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'OUTPUT_PARTICLES_FMT') THEN
              
              this%output_particles_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'OUTPUT_PARTICLES_FREQ_STEP') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200) &
                   this%output_particles_freq_step
              
           ELSE IF (carg == 'OUTPUT_PARTICLES_FREQ_TIME') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200) &
                   this%output_particles_freq_time
              
           ELSE  IF (carg == 'OUTPUT_CONFORMATION_FILE') THEN
              
              this%output_conformation_file = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'OUTPUT_CONFORMATION_FMT') THEN
              
              this%output_conformation_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'OUTPUT_CONFORMATION_FREQ_STEP') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200) &
                   this%output_conformation_freq_step
              
           ELSE IF (carg == 'OUTPUT_CONFORMATION_FREQ_TIME') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200) &
                   this%output_conformation_freq_time
              
           ELSE  IF (carg == 'COLLOID_FILE') THEN
              
              this%colloid_file = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'COLLOID_FMT') THEN
              
              this%colloid_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'COLLOID_FREQ_STEP') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200)  this%colloid_freq_step
              
           ELSE IF (carg == 'COLLOID_FREQ_TIME') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200)  this%colloid_freq_time
     
           ELSE  IF (carg == 'STATISTIC_RELAX_FILE') THEN
              
              this%statistic_relax_file = TRIM(ADJUSTL(cvalue))
             
           ELSE IF (carg == 'STATISTIC_RELAX_FMT') THEN
              
              this%statistic_relax_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'STATISTIC_RELAX_FREQ_STEP') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200)  &
                   this%statistic_relax_freq_step
              
           ELSE  IF (carg == 'STATISTIC_FILE') THEN
              
              this%statistic_file = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'STATISTIC_FMT') THEN
              
              this%statistic_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'STATISTIC_FREQ_STEP') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200)  this%statistic_freq_step
              
           ELSE IF (carg == 'STATISTIC_FREQ_TIME') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200)  this%statistic_freq_time
             
           ELSE  IF (carg == 'BOUNDARY_FILE') THEN
              
              this%boundary_file = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'BOUNDARY_FMT') THEN
              
              this%boundary_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'BOUNDARY_FREQ_STEP') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200)  this%boundary_freq_step
      
           ELSE IF (carg == 'BOUNDARY_FREQ_TIME') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200)  this%boundary_freq_time
              
           ELSE IF (carg == 'RESTART_PHYSICS_FILE') THEN
              
              this%restart_physics_file = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'RESTART_PHYSICS_FMT') THEN
              
              this%restart_physics_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE  IF (carg == 'RESTART_PARTICLES_RELAX_FILE') THEN
              
              this%restart_particles_relax_file = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'RESTART_PARTICLES_RELAX_FMT') THEN
              
              this%restart_particles_relax_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE  IF (carg == 'RESTART_PARTICLES_FILE') THEN
              
              this%restart_particles_file = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'RESTART_PARTICLES_FMT') THEN
              
              this%restart_particles_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE  IF (carg == 'RESTART_CONFORMATION_FILE') THEN
              
              this%restart_conformation_file = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'RESTART_CONFORMATION_FMT') THEN
              
              this%restart_conformation_fmt = TRIM(ADJUSTL(cvalue))
              
           ELSE IF (carg == 'RESTART_FREQ_STEP') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200)  this%restart_freq_step
              
           ELSE IF (carg == 'RESTART_FREQ_TIME') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200)  this%restart_freq_time
              
           ELSE IF (carg == 'RESTART_FREQ_TIME_WALL') THEN
              
              READ(cvalue,*,IOSTAT=ios, ERR=200)  this%restart_freq_time_wall
              
           END IF
           
        END DO
        
        !----------------------------------------------------
        ! End of file
        !----------------------------------------------------
        
       
        !----------------------------------------------------
        ! Something went wrong
        !----------------------------------------------------
        
200     CONTINUE
        
        WRITE(cbuf,'(A,I5,2A)') 'Error reading line: ',iline,&
             &  ' of file: ',this%io_config_file(1:ilenio)
	
        ilen = LEN_TRIM(cbuf)
        PRINT *,'io_read_io_config : ',cbuf(1:ilen)
        
        stat_info = -1
        
        
9999    CONTINUE
	
	!----------------------------------------------------
        ! Close file
      	!----------------------------------------------------
	
        CLOSE(this%io_config_unit)
        
#ifdef __DEBUG
        IF ( debug_threshold > 1) THEN
           PRINT *, "io_read_io_config : ended"
        END IF
#endif
        
        RETURN        
        
      END SUBROUTINE io_read_io_config
      
