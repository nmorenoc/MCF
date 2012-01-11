      SUBROUTINE io_write_conformation(this,rank,step,&
           d_particles, num_part,stat_info)
        !-------------------------------------------------------------
        ! Subroutine  :  io_write_conformation
        !-------------------------------------------------------------
        !
        ! Purpose     :  Writing particles' conformation into files,
        !                in case of Non-Newtonian fluids,
        !                taking a Particles object as an argument.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   :V0.1 25.08.2010, If the eigen_dynamics variable
        !              is .TRUE. the eigenvalues and 
        !              eigenvectors are written to the output file.
        !              If it is .FALSE. the conformation tensor is 
        !              written to the output file. (Adolfo)
        !              
        !              V0.1 22.04.2010, pressure tensor now is 
        !              available to output. (Xin)
        !          
        !              V0.1 03.08.2009, original version.
        !
        !-------------------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !-------------------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !
        ! this           : an object of Marching Class.
        ! rank           : rank of process.
        ! step           : index of current step.
        ! d_particles    : an object of Particles Class
        ! num_part       : first num_part needed to be written.
        ! stat_info      : return flag of status.
        !----------------------------------------------------
        
        !----------------------------------------------------
      	! Modules
      	!----------------------------------------------------
        
        USE ppm_module_user_io
        
        
        !----------------------------------------------------
        ! Arguments     
      	!----------------------------------------------------
        
        TYPE(IO),INTENT(IN)                     :: this
        INTEGER, INTENT(IN)                     :: rank
        INTEGER, INTENT(IN)                     :: step
        TYPE(Particles), INTENT(IN)             :: d_particles
        INTEGER, INTENT(IN)                     :: num_part        
        INTEGER,INTENT(OUT)                     :: stat_info
        
      	!----------------------------------------------------
      	! Local variables 
      	!----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        
        LOGICAL                                 :: read_external
        LOGICAL                                 :: eigen_dynamics
        REAL(MK), DIMENSION(:,:)  , POINTER     :: x
        INTEGER , DIMENSION(:)    , POINTER     :: pid
        REAL(MK), DIMENSION(:,:)  , POINTER     :: vgt
        REAL(MK), DIMENSION(:,:,:), POINTER     :: pt
        REAL(MK), DIMENSION(:,:)  , POINTER     :: ct
        REAL(MK), DIMENSION(:,:)  , POINTER     :: eval
        REAL(MK), DIMENSION(:,:)  , POINTER     :: evec
        REAL(MK), DIMENSION(:,:)  , POINTER     :: output

        CHARACTER(LEN=MAX_CHAR)                 :: file_name
        INTEGER                                 :: file_fmt
        INTEGER                                 :: file_unit
        
        INTEGER                                 :: dim
        INTEGER                                 :: data_dim
        INTEGER                                 :: cur_dim
        INTEGER                                 :: i, j
        
        CHARACTER(LEN=MAX_CHAR)                 :: cbuf
        INTEGER					:: clen
        
	!----------------------------------------------------
      	! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(x)
        NULLIFY(pid)
        NULLIFY(vgt)
        NULLIFY(pt)
        NULLIFY(ct)
        NULLIFY(eval)
        NULLIFY(evec)
        NULLIFY(output)
        
        read_external  = &
             control_get_read_external(this%ctrl,stat_info_sub)
        dim            = &
             physics_get_num_dim(this%phys,stat_info_sub)
        eigen_dynamics = &
             physics_get_eigen_dynamics(this%phys,stat_info_sub)
        
        CALL particles_get_x(d_particles,x,num_part,stat_info_sub)
        CALL particles_get_pid(d_particles,pid,num_part,stat_info_sub)
        CALL particles_get_vgt(d_particles,vgt,num_part,stat_info_sub)
        CALL particles_get_pt(d_particles,pt,num_part,stat_info_sub)
        IF (.NOT.(eigen_dynamics)) THEN
           CALL particles_get_ct(d_particles,ct,num_part,stat_info_sub)
        ELSE
           CALL particles_get_eval(d_particles,eval,num_part,stat_info_sub)
           CALL particles_get_evec(d_particles,evec,num_part,stat_info_sub)           
        ENDIF

        IF (.NOT.(eigen_dynamics)) THEN
           data_dim = dim + 1 + dim**2 + dim**2 + dim**2
        ELSE
           data_dim = dim + 1 + dim**2 + dim**2 + dim + dim**2
        ENDIF
        
        ALLOCATE(output(data_dim,num_part))
        
        output(1:dim,1:num_part) = x(1:dim,1:num_part)
        
        cur_dim = dim + 1 
        
        output(cur_dim,1:num_part) = pid(1:num_part)
        
        cur_dim = cur_dim+ 1 
        
        output(cur_dim:cur_dim+dim**2-1,1:num_part) = &
             vgt(1:dim**2,1:num_part)

        cur_dim = cur_dim + dim**2 

        IF (.NOT.(eigen_dynamics)) THEN
           output(cur_dim:cur_dim+dim**2-1,1:num_part) = &
                ct(1:dim**2,1:num_part)
           
           cur_dim = cur_dim + dim**2 
        ELSE
           output(cur_dim:cur_dim+dim-1,1:num_part) = &
                eval(1:dim,1:num_part)

           cur_dim = cur_dim + dim 

           output(cur_dim:cur_dim+dim**2-1,1:num_part) = &
                evec(1:dim**2,1:num_part)
           
           cur_dim = cur_dim + dim**2 
        ENDIF

        do i = 1, dim
           do j = 1, dim
              output(cur_dim + j + dim*(i-1) - 1,1:num_part) = &
                   pt(i,j,1:num_part)
           enddo
        enddo

        cur_dim = cur_dim + dim**2 

        !----------------------------------------------------
      	! Define the output file name for this time step.
        ! If we have read particles from file,
        ! the output should include "_r" to distingusih.
      	!----------------------------------------------------
        
        IF(read_external) THEN
           WRITE(file_name,'(2A,I8.8,A)') &
                TRIM(this%output_conformation_file),'_r',step,'.out'
        ELSE
           WRITE(file_name,'(A,I8.8,A)') &
                TRIM(this%output_conformation_file),step,'.out'
        END IF
        
      	!----------------------------------------------------
      	! Define format of output file.
      	!----------------------------------------------------
        
        IF (this%output_conformation_fmt(1:1) .EQ. 'f' .OR. &
             this%output_conformation_fmt(1:1) .EQ. 'F') THEN
           file_fmt = ppm_param_io_ascii
        ELSE
           file_fmt = ppm_param_io_binary
        END IF
        
        file_unit = this%output_conformation_unit        
        
        !----------------------------------------------------
        ! Open ppm I/O unit for centralized I/O.
        !----------------------------------------------------
        
        CALL ppm_io_open(file_unit,file_name,&
             ppm_param_io_write, ppm_param_io_replace, &
             file_fmt,ppm_param_io_centralized,stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN          
           PRINT *, 'io_write_conformation : ', &
                'Failed to open unit !'
           stat_info = -1
           GOTO 9999
        END IF
        
        WRITE(cbuf,'(A1,I2,A6)') '(', data_dim ,'E16.8)'
        clen = LEN_TRIM(cbuf)         
        
        CALL ppm_io(file_unit,output,ppm_param_io_write,&
             ppm_param_io_concat,IOFMT=cbuf(1:clen), STAT=stat_info_sub)
        
        IF (stat_info_sub /= 0) THEN          
           PRINT *, 'io_write_conformation : ',&
                "Writing conformation failed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Close file.
        !----------------------------------------------------
        
        CALL ppm_io_close(file_unit,stat_info_sub)
        
        !----------------------------------------------------
        ! Print out confirm.
        !----------------------------------------------------
        
        IF (rank == 0) THEN
           IF (.NOT.(eigen_dynamics)) THEN
              WRITE(cbuf,'(2A)') 'Conformation tensor written to ',&
                   TRIM(file_name)
           ELSE
              WRITE(cbuf,'(2A)') 'Eigenvectors and eigenvalues written to ',&
                   TRIM(file_name)
           ENDIF

           PRINT *, TRIM(cbuf)
        END IF
        
        !----------------------------------------------------
        ! Return.
        !----------------------------------------------------
        
9999    CONTINUE
        
        IF (ASSOCIATED(x)) THEN
           DEALLOCATE(x)
        END IF
        
        IF (ASSOCIATED(pid)) THEN
           DEALLOCATE(pid)
        END IF
        
        IF (ASSOCIATED(vgt)) THEN
           DEALLOCATE(vgt)
        END IF
        
        IF (ASSOCIATED(pt)) THEN
           DEALLOCATE(pt)
        END IF
        
        IF (ASSOCIATED(ct)) THEN
           DEALLOCATE(ct)
        END IF
        
        IF (ASSOCIATED(output)) THEN
           DEALLOCATE(output)
        END IF
        
        RETURN	
        
      END SUBROUTINE io_write_conformation
