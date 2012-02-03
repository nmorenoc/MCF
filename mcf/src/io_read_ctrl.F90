      SUBROUTINE io_read_ctrl(this,ctrl,stat_info)
        !----------------------------------------------------
        ! Subroutine  :  io_read_ctrl
        !----------------------------------------------------
        !
        ! Purpose     : Reading control parameters from
        !               control file.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 03.03.2009, original version.
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
        ! Arguments.
        !----------------------------------------------------
        
        TYPE(IO), INTENT(INOUT)         :: this
        TYPE(Control), INTENT(INOUT)    :: ctrl
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
      	! Local variables.
      	!----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                        	:: i,j,idx
        INTEGER                        	:: ilen,ios
        INTEGER                         :: iline,ilenctrl
        CHARACTER(LEN=MAX_CHAR)         :: cbuf,carg,cvalue
        LOGICAL                        	:: lExist  
        
        INTEGER                         :: debug_flag
        LOGICAL                         :: relax_run
        LOGICAL                         :: read_external
        INTEGER                         :: kernel_type
        LOGICAL                         :: symmetry
        INTEGER                         :: rhs_density_type
        LOGICAL                         :: dynamic_density_ref
        INTEGER                         :: stateEquation_type
        LOGICAL                         :: Newtonian
        LOGICAL                         :: Brownian
        INTEGER                         :: random_seed
        INTEGER                         :: rhs_force_type
        LOGICAL                         :: pp_interact_cc
        LOGICAL                         :: pp_interact_cw
        INTEGER                         :: cc_lub_type
        INTEGER                         :: cc_repul_type
        INTEGER                         :: cw_lub_type
        INTEGER                         :: cw_repul_type
        LOGICAL                         :: p_energy
        LOGICAL                         :: flow_v_fixed
        INTEGER                         :: integrate_type
        INTEGER                         :: adaptive_dt
        INTEGER                         :: write_output
        INTEGER                         :: write_restart

#ifdef __DEBUG
        INTEGER                         :: debug_threshold
#endif

	!----------------------------------------------------
      	! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info    = 0
        stat_info_sub = 0
        
        cbuf   = ""
        carg   = ""
        cvalue = ""
        
#ifdef __DEBUG
        
        debug_threshold = 0
        
        IF ( debug_threshold > 1) THEN           
           PRINT *, "io_read_ctrl : starting "          
        END IF
#endif
        
        
        !----------------------------------------------------
        ! Check if the name of ctrl file if empty.
      	!----------------------------------------------------
       
        ilenctrl = LEN_TRIM(this%ctrl_file)
        
        IF (ilenctrl < 1) THEN
           PRINT *,'io_read_ctrl : ',&
                'No Ctrl file name given!'
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check if the ctrl file exists.
        !----------------------------------------------------
        
        INQUIRE(FILE=this%ctrl_file,EXIST=lExist)
        IF ( .NOT.lExist ) THEN
           WRITE(cbuf,'(2A)')'No such file: ', &
                this%ctrl_file(1:ilenctrl)
           PRINT *, 'io_read_ctrl : ', cbuf
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Open control file for reading.
      	!----------------------------------------------------
        
        OPEN(this%ctrl_unit,FILE=this%ctrl_file, &
             IOSTAT=ios,ACTION='READ')
        
        IF (ios /= 0) THEN 
           WRITE(cbuf,'(2A)')'Failed to open file: ',&
                this%ctrl_file(1:ilenctrl)
           PRINT *,'io_read_ctrl : ', cbuf
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Scan the file line by line.
        !----------------------------------------------------
        
        iline = 0
	
      	DO 
          
          !--------------------------------------------------
          ! Increase line counter.
          !--------------------------------------------------
          
          iline = iline + 1
          
          !--------------------------------------------------
          ! Read current line.
          !--------------------------------------------------
          
          READ(this%ctrl_unit,'(A)',END=9999,ERR=200) cbuf
          
          !--------------------------------------------------
          ! Count length of current string read.
          !--------------------------------------------------
          
          ilen = LEN_TRIM(cbuf)
          
          !--------------------------------------------------
          ! Extensive Debug: print each line of Ctrl file.
          !--------------------------------------------------
          
#ifdef __DEBUG
          
          IF ( debug_threshold > 3 ) THEN             
             PRINT *,'io_read_ctrl', cbuf(1:ilen)
          END IF
#endif          
          
          !--------------------------------------------------
          ! Skip empty line or 
          ! comment line with symbol "#" at first place.
          !--------------------------------------------------
          
          IF (ilen < 1 .OR. cbuf(1:1) == '#' ) THEN
             CYCLE
          END IF
          
          !--------------------------------------------------
          ! Remove spaces of the line being read.
          !--------------------------------------------------
          
          j = 0
          DO i=1,ilen
             IF (cbuf(i:i) /=' ' .AND. cbuf(i:i) /= '\t' ) THEN
                j = j + 1
                cbuf(j:j) = cbuf(i:i)
             END IF
          END DO
          
          !--------------------------------------------------
          ! Update length of string.
          !--------------------------------------------------
          
          ilen = j
          
          !--------------------------------------------------
          ! After removing spaces,
          ! skip comment line with symbol "#" at first place.
          !--------------------------------------------------
          
          IF (cbuf(1:1) == '#' ) THEN
             CYCLE
          END IF
          
          !--------------------------------------------------
          ! Find position of "=".
          !--------------------------------------------------
          
          idx = INDEX(cbuf,'=')
          
          !--------------------------------------------------
          ! Exit if "=" is missing.
          !--------------------------------------------------
          
          IF (idx < 0) THEN
             WRITE(cbuf,'(A,I5)')'Incorrect line without = : ',iline
             Print *,'io_read_ctrl : ', cbuf
             stat_info = -1
             GOTO 9999
          END IF
          
          !--------------------------------------------------
          ! Get name of argument and its value.
          !--------------------------------------------------
          
          carg   = ADJUSTL(cbuf(1:idx-1))
          cvalue = ADJUSTL(cbuf(idx+1:ilen))
          
          !--------------------------------------------------
          ! Convert name of argument to upper case.
          !--------------------------------------------------
          
          CALL tool_uppercase(this%io_tool,carg,idx-1,stat_info)
          
#ifdef __DEBUG
          
          IF ( debug_threshold > 2) THEN
             PRINT *, 'io_read_ctrl : ', carg
             PRINT *, 'io_read_ctrl : ', cvalue
          END IF
#endif     
          
          IF (carg == 'JOB_NAME') THEN
             
             CALL control_set_job_name(ctrl,&
                  cvalue(1:LEN_TRIM(cvalue)),stat_info_sub)
             
          ELSE IF (carg == 'JOB_SUBMIT_DATE') THEN
             
             CALL control_set_job_submit_date(ctrl,&
                  cvalue(1:LEN_TRIM(cvalue)),stat_info_sub)
             
          ELSE IF (carg == 'MCF_PHYSICS_CONFIG_FILE') THEN
             
             this%physics_config_file = cvalue
             
          ELSE IF (carg == 'MCF_IO_CONFIG_FILE') THEN
             
             this%io_config_file = cvalue             
             
          ELSE IF (carg == 'DEBUG_FLAG') THEN
	     
             READ(cvalue,*,IOSTAT=ios,ERR=200) debug_flag
             CALL control_set_debug_flag(ctrl,&
                  debug_flag,stat_info_sub)
             
          ELSE IF (carg == 'RELAX_RUN') THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) relax_run
             CALL control_set_relax_run(ctrl,&
                  relax_run,stat_info_sub)
             
          ELSE IF (carg == 'READ_EXTERNAL') THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) read_external
             CALL control_set_read_external(ctrl,&
                  read_external,stat_info_sub)
             
          ELSE IF (carg == 'KERNEL_TYPE') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) kernel_type
             CALL control_set_kernel_type(ctrl,&
                  kernel_type,stat_info_sub)
             
          ELSE IF (carg == 'SYMMETRY') THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) symmetry
             CALL control_set_symmetry(ctrl,symmetry,stat_info_sub)
             
          ELSE IF (carg == 'RHS_DENSITY_TYPE') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) rhs_density_type
             CALL  control_set_rhs_density_type(ctrl,&
                  rhs_density_type,stat_info_sub)
             
          ELSE IF (carg == 'DYNAMIC_DENSITY_REF') THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) dynamic_density_ref
             CALL control_set_dynamic_density_ref(ctrl,&
                  dynamic_density_ref,stat_info_sub)
             
          ELSE IF (carg == 'STATEEQUATION_TYPE') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) stateEquation_type
             CALL  control_set_stateEquation_type(ctrl,&
                  stateEquation_type,stat_info_sub)
             
          ELSE IF (carg == 'NEWTONIAN') THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) Newtonian
             CALL control_set_Newtonian(ctrl,Newtonian,stat_info_sub)
             
          ELSE IF (carg == 'BROWNIAN') THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) Brownian
             CALL control_set_Brownian(ctrl,Brownian,stat_info_sub)
             
          ELSE IF (carg == 'RANDOM_SEED') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) random_seed
             CALL  control_set_random_seed(ctrl,&
                  random_seed,stat_info_sub)
             
          ELSE IF (carg == 'RHS_FORCE_TYPE') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) rhs_force_type
             CALL  control_set_rhs_force_type(ctrl,&
                  rhs_force_type,stat_info_sub)
             
          ELSE IF (carg == 'PP_INTERACT_CW') THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) pp_interact_cw
             CALL  control_set_pp_interact_cw(ctrl,&
                  pp_interact_cw,stat_info_sub)
             
          ELSE IF (carg == 'PP_INTERACT_CC') THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) pp_interact_cc
             CALL  control_set_pp_interact_cc(ctrl,&
                  pp_interact_cc,stat_info_sub)
             
          ELSE IF (carg == 'CC_LUB_TYPE') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cc_lub_type
             CALL  control_set_cc_lub_type(ctrl,&
                  cc_lub_type,stat_info_sub)
             
          ELSE IF (carg == 'CC_REPUL_TYPE') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cc_repul_type
             CALL  control_set_cc_repul_type(ctrl,&
                  cc_repul_type,stat_info_sub)
             
          ELSE IF (carg == 'CW_LUB_TYPE') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cw_lub_type
             CALL  control_set_cw_lub_type(ctrl,&
                  cw_lub_type,stat_info_sub)
             
          ELSE IF (carg == 'CW_REPUL_TYPE') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) cw_repul_type
             CALL  control_set_cw_repul_type(ctrl,&
                  cw_repul_type,stat_info_sub)
             
          ELSE IF (carg == 'P_ENERGY') THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) p_energy
             CALL control_set_p_energy(ctrl,&
                  p_energy,stat_info_sub)
     
          ELSE IF (carg == 'FLOW_V_FIXED') THEN
             
             READ(cvalue,'(L)',IOSTAT=ios,ERR=200) flow_v_fixed
             CALL control_set_flow_v_fixed(ctrl,&
                  flow_v_fixed,stat_info_sub)
             
          ELSE IF (carg == 'INTEGRATE_TYPE') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) integrate_type
             CALL control_set_integrate_type(ctrl,&
                  integrate_type,stat_info_sub)
             
          ELSE IF (carg == 'ADAPTIVE_DT') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) adaptive_dt
             CALL control_set_adaptive_dt(ctrl,&
                  adaptive_dt,stat_info_sub)

          ELSE IF (carg == 'WRITE_OUTPUT') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) write_output
             CALL control_set_write_output(ctrl,&
                  write_output,stat_info_sub)
             
          ELSE IF (carg == 'WRITE_RESTART') THEN
             
             READ(cvalue,*,IOSTAT=ios,ERR=200) write_restart
             CALL control_set_write_restart(ctrl,&
                  write_restart,stat_info_sub)
          END IF
          
        END DO
       
        !----------------------------------------------------
        ! End of file.
      	!----------------------------------------------------
        
        
      	!----------------------------------------------------
      	! Something went wrong.
      	!----------------------------------------------------
        
200     CONTINUE
        
        WRITE(cbuf,'(A,I5,2A)') &
             'Error reading line: ',iline,&
             ' of file: ',this%ctrl_file(1:ilenctrl)
	
        ilen = LEN_TRIM(cbuf)
        PRINT *,'io_read_ctrl : ',cbuf(1:ilen)
        
        stat_info = -1
        
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Close file.
        !----------------------------------------------------
        CLOSE(this%ctrl_unit)

#ifdef __DEBUG
        
        IF ( debug_threshold > 1) THEN
           PRINT *, "io_read_ctrl : ended "
	END IF
 
#endif

        RETURN        
 
      END SUBROUTINE io_read_ctrl
