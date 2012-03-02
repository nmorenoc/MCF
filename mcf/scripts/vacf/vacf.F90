      PROGRAM vacf
        !----------------------------------------------------
        ! Calculate Velocity AutoCorrelation Function of 
        ! a colloid.
        !
        ! Input 
        !
        ! num_dim  : number of dimension
        ! num_line : number of lines in the file.
        ! 
        ! Output   :
        !             dt, msd
        !----------------------------------------------------
     
        IMPLICIT NONE
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------        
        
        INTEGER, PARAMETER              :: MK= KIND(1.0D0)
        INTEGER                         :: stat_info
        
        INTEGER                         :: n_arg
        CHARACTER(256)                  :: cbuf
        INTEGER                         :: num_dim
        INTEGER                         :: num_line
        INTEGER                         :: num_dt
        REAL(MK), DIMENSION(:),POINTER  :: time
        REAL(MK), DIMENSION(3)          :: r
        REAL(MK), DIMENSION(:,:),POINTER:: v
        
        
        CHARACTER(256)                  :: filename
        INTEGER                         :: fileunit
        INTEGER                         :: flength
        
        LOGICAL                         :: lExist
        LOGICAL                         :: lOpened
        
        INTEGER                         :: i
        
        INTEGER                         :: n_dt,n0,n1
        REAL(MK), DIMENSION(:), POINTER :: vacf_c
        REAL(MK), DIMENSION(:), POINTER :: D
        REAL(MK), DIMENSION(:),POINTER  :: dt
        
        REAL(MK)                        :: ddt        
        CHARACTER(256)                  :: vacfname
        INTEGER                         :: vacfunit
        INTEGER                         :: num_thread
        INTEGER                         :: thread_num
        INTEGER                         :: OMP_GET_THREAD_NUM
        INTEGER                         :: OMP_GET_NUM_THREADS
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        NULLIFY(time)
        NULLIFY(v)
        NULLIFY(vacf_c)
        NULLIFY(D)
        NULLIFY(dt)

        !----------------------------------------------------
        ! Get input parameters from command line.
        !----------------------------------------------------
     
        n_arg = IARGC()
        
        IF (n_arg > 1 ) THEN
           CALL GETARG(1,cbuf)
           READ(cbuf,*), num_dim
           CALL GETARG(2,cbuf)
           READ(cbuf,*), num_line
        ELSE
           PRINT *, "command line argument not enough !"
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check input paramters.
        !----------------------------------------------------
  
        PRINT *, "num_dim : ", num_dim
        PRINT *, "num_line : ", num_line
        
        IF ( num_dim < 2 .OR. num_dim > 3 ) THEN
           PRINT *, "num_dim not 2 or 3 !"
           GOTO 9999
        END IF
        
        IF ( num_line < 0 ) THEN
           PRINT *, "num_line < 0 !"
           GOTO 9999
        END IF
        
        
        !------------------------------------
        ! Check if the particle file exists. 
        !------------------------------------
        
        filename = TRIM("mcf_colloid.dat")
        flength  = LEN_TRIM(filename)
        
        lExist = .TRUE.
        INQUIRE(FILE=filename,EXIST=lExist)
        
        IF (.NOT.lExist) THEN
           WRITE(cbuf,'(2A)')'No such file: ', &
                filename(1:flength)
           PRINT *, cbuf
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Open the input file.
        !----------------------------------------------------
        
        stat_info = 0
        OPEN(fileunit,FILE=filename,&
             IOSTAT=stat_info,ACTION='READ')
        
        IF (stat_info /= 0) THEN
           WRITE(cbuf,'(2A)')'Failed to open file: ',&
                filename(1:flength)
           PRINT *, cbuf
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Allocate memory to record input data.
        !----------------------------------------------------
        
        ALLOCATE(time(0:num_line))
        ALLOCATE(v(num_dim,0:num_line))
        ALLOCATE(dt(0:num_line))
        
        !----------------------------------------------------
        ! Scan the file and save data.
        !----------------------------------------------------
        
        DO i = 0, num_line
           
           READ(fileunit,*,END=9999,ERR=200) &
                time(i), r(1:num_dim), v(1:num_dim,i)
           
        END DO
        
        !----------------------------------------------------
        ! Truncate the number of dt, since the biggest dt
        !  won't have enough statistics.
        !----------------------------------------------------
        
        num_dt = num_line * 4 / 5
        
        ALLOCATE(vacf_c(0:num_dt))
        ALLOCATE(D(-1:num_dt))
        
        D(-1)  = 0.0_MK
        D(0)   = 0.0_MK
        
        !----------------------------------------------------
        ! Use OpenMP to calcuate the MSD.
        !----------------------------------------------------
        
        !$OMP PARALLEL PRIVATE(thread_num)
        
        num_thread = OMP_GET_NUM_THREADS()
        thread_num = OMP_GET_THREAD_NUM()
        
        PRINT *, "num_thread, thread_num : ", &
             num_thread, thread_num
        
        !$OMP DO
        
        DO n_dt = 0, num_dt
           
           !-------------------------------------------------
           ! Reset vacf and record dt.
           !-------------------------------------------------
           
           vacf_c(n_dt) = 0.0_MK
           dt(n_dt) = time(n_dt) - time(0)
           
           !-------------------------------------------------
           ! r data starts from line 0.
           !-------------------------------------------------
           
           DO n0 = 0, num_line-n_dt
              
              !----------------------------------------------
              ! reference v0 is n0, v1 is at n0+n_dt,
              ! calculate v(n0) dot product v(n1).
              !----------------------------------------------
              
              n1 = n0 + n_dt
              
              vacf_c(n_dt) = vacf_c(n_dt) + &
                   v(1,n0) * v(1,n1) + &
                   v(2,n0) * v(2,n1)
              
              IF ( num_dim == 3 ) THEN
                 
                 vacf_c(n_dt) = vacf_c(n_dt) + &
                      v(3,n0) * v(3,n1)
                 
              END IF
              
              
           END DO ! n0 = 0, num_line-n_dt
           
           !-------------------------------------------------
           ! Calcuate average vacf at dt = n_dt.          
           !-------------------------------------------------
           
           vacf_c(n_dt) = vacf_c(n_dt) / REAL((num_line-n_dt+1),MK)
           
        END DO ! n_dt = 0, num_dt

        !$OMP END DO

        !$OMP END PARALLEL

        !----------------------------------------------------
        ! Calcualte diffusion coefficent.
        !----------------------------------------------------
        
        ddt = time(1) - time(0)

        DO n_dt = 1, num_dt
           D(n_dt) = D(n_dt-1) + &
                (vacf_c(n_dt-1) + vacf_c(n_dt)) * ddt /2.0
        END DO
        
        
        !----------------------------------------------------
        ! output dt, vacf and diffusion coefficient into file.
        !----------------------------------------------------
        
        vacfname =TRIM("vacf.dat")
        
        OPEN(UNIT=vacfunit,FILE=vacfname,&
             STATUS="REPLACE",FORM="FORMATTED",ACTION="WRITE", &
             IOSTAT=stat_info)
        IF (stat_info /= 0) THEN
           PRINT *, "Writting vacf file has problem !"
        END IF
        
        DO i = 0, num_dt
           
           WRITE(vacfunit,FMT='(3(E16.8))',IOSTAT=stat_info) &
                dt(i), vacf_c(i), D(i)/num_dim
           
        END DO
        
        PRINT *, "VACF written to ", TRIM(vacfname)
        GOTO 9999
        
200     CONTINUE
        
        PRINT *, "reading line ", i, " has error !"
        
9999    CONTINUE
        
        INQUIRE(UNIT=fileunit,OPENED=lOpened)
        IF (lopened) THEN
           CLOSE(fileunit)
        END IF
        
        INQUIRE(UNIT=vacfunit,OPENED=lOpened)
        IF (lopened) THEN
           CLOSE(vacfunit)
        END IF
        
        IF(ASSOCIATED(time)) THEN
           DEALLOCATE(time)
        END IF
        
        IF(ASSOCIATED(v)) THEN
           DEALLOCATE(v)
        END IF
        
        IF(ASSOCIATED(vacf_c)) THEN
           DEALLOCATE(vacf_c)
        END IF
        
        IF(ASSOCIATED(D)) THEN
           DEALLOCATE(D)
        END IF
        
        IF(ASSOCIATED(dt)) THEN
           DEALLOCATE(dt)
        END IF
        
      END PROGRAM vacf
      
      
      
