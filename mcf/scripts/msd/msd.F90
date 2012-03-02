      PROGRAM msd
        !----------------------------------------------------
        ! Calculate Mean Square Displacement of a colloid
        ! (2D/3D).
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
        REAL(MK),DIMENSION(:), POINTER  :: time
        REAL(MK),DIMENSION(:,:),POINTER :: r
        REAL(MK),DIMENSION(3)           :: v
        
        CHARACTER(256)                  :: filename
        INTEGER                         :: fileunit
        INTEGER                         :: flength
        
        LOGICAL                         :: lExist
        LOGICAL                         :: lOpened
        
        INTEGER                         :: i
        
        INTEGER                         :: n_dt,n0,n1
        REAL(MK), DIMENSION(:), POINTER :: msd_c
        REAL(MK), DIMENSION(3)          :: rij
        REAL(MK),DIMENSION(:), POINTER  :: dt
        
        CHARACTER(256)                  :: msdname
        INTEGER                         :: msdunit
        INTEGER                         :: thread_num
        INTEGER                         :: num_thread
        INTEGER                         :: OMP_GET_THREAD_NUM
        INTEGER                         :: OMP_GET_NUM_THREADS
    
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        
        NULLIFY(time)
        NULLIFY(r)
        NULLIFY(msd_c)
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
        
        !----------------------------------------------------
        ! Check if the given file exists. 
        !----------------------------------------------------
        
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
        OPEN(fileunit,FILE=filename,IOSTAT=stat_info, &
             ACTION='READ')
        
        IF ( stat_info /= 0 ) THEN
           WRITE(cbuf,'(2A)')'Failed to open file: ', &
                filename(1:flength)
           PRINT *, cbuf
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Allocate memory to record input data.
        !----------------------------------------------------
        
        ALLOCATE(time(0:num_line))
        ALLOCATE(r(num_dim,0:num_line))
        
        !----------------------------------------------------
        ! Scan the file and save data.
        !----------------------------------------------------
        
        DO i = 0, num_line
           
           READ(fileunit,*,END=9999,ERR=200) &
                time(i), r(1:num_dim,i), v(1:num_dim)
           
        END DO ! i = 0, num_line
        
        !----------------------------------------------------
        ! Truncate the number of dt, since the biggest dt
        ! won't have enough statistics.
        !----------------------------------------------------
        
        num_dt = num_line * 4 / 5
        
        ALLOCATE(dt(0:num_dt))
        ALLOCATE(msd_c(0:num_dt))
        
        dt(0)    = 0.0_MK
        msd_c(0) = 0.0_MK
        
        
        !----------------------------------------------------
        ! Use OpenMP to calcualte the MSD in parallel.
        !----------------------------------------------------
        
        !$OMP PARALLEL
        
        !num_thread = OMP_GET_NUM_THREADS()
        !thread_num = OMP_GET_THREAD_NUM()
        
        PRINT *, "num_thread, thread_num : ", &
             num_thread, thread_num
        
        
        !$OMP DO
        
        !----------------------------------------------------
        ! Loop over different sizes of dt, from 1 to num_dt.
        !----------------------------------------------------
        
        DO n_dt = 1, num_dt
           
           !-------------------------------------------------
           ! Reset msd and record dt.
           !-------------------------------------------------
           
           msd_c(n_dt) = 0.0_MK
           dt(n_dt)    = time(n_dt) - time(0)
           
           !-------------------------------------------------
           ! r data starts from 0 line.
           !-------------------------------------------------
           
           DO n0 = 0, num_line - n_dt
              
              !----------------------------------------------
              ! reference r0 is n0, r1 is at n0+n_dt,
              ! calculate |r1-r0|**2 at n1.
              !----------------------------------------------
              
              n1 = n0 + n_dt
              
              rij(1:num_dim) = r(1:num_dim,n1) - &
                   r(1:num_dim,n0)
              
              msd_c(n_dt) = msd_c(n_dt) + &
                   rij(1)**2 + rij(2)**2

              IF ( num_dim == 3 ) THEN
                 msd_c(n_dt) = msd_c(n_dt) + &
                      rij(3)**2
              END IF
              
           END DO ! n0 = 0, num_line - n_dt
           
           !-------------------------------------------------
           ! Calcuate average msd at dt = n_dt.          
           !-------------------------------------------------
           
           msd_c(n_dt) = msd_c(n_dt) / REAL((num_line-n_dt+1),MK)
           
        END DO ! n_dt = 1, num_dt
        
        !$OMP END DO
        
        !$OMP END PARALLEL
        
        !----------------------------------------------------
        ! output dt and msd into file.
        !----------------------------------------------------
        
        msdname =TRIM("msd.dat")
        
        OPEN(UNIT=msdunit,FILE=msdname,&
             STATUS="REPLACE",FORM="FORMATTED",ACTION="WRITE", &
             IOSTAT=stat_info)
        
        IF ( stat_info /= 0 ) THEN
           PRINT *, "Writting msd file has problem !"
           GOTO 9999
        END IF
        
        
        DO i = 0, num_dt
           
           WRITE(msdunit,FMT='(2(E16.8))',IOSTAT=stat_info) &
                dt(i), msd_c(i)
           
        END DO
        
        PRINT *, "MSD written to ", TRIM(msdname)
        
        GOTO 9999
        
200     CONTINUE
        
        PRINT *, "reading line ", i, " has error !"
        
9999    CONTINUE
        
        INQUIRE(UNIT=fileunit,OPENED=lOpened)
        IF (lopened) THEN
           CLOSE(fileunit)
        END IF
        
        INQUIRE(UNIT=msdunit,OPENED=lOpened)
        IF (lopened) THEN
           CLOSE(msdunit)
        END IF
        
        IF(ASSOCIATED(time)) THEN
           DEALLOCATE(time)
        END IF
        
        IF(ASSOCIATED(r)) THEN
           DEALLOCATE(r)
        END IF
        
        IF(ASSOCIATED(msd_c)) THEN
           DEALLOCATE(msd_c)
        END IF
        
        IF(ASSOCIATED(dt)) THEN
           DEALLOCATE(dt)
        END IF
        
      END PROGRAM msd
      
      
      
