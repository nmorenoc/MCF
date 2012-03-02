      PROGRAM unfold_r
        
        IMPLICIT NONE
        
        INTEGER, PARAMETER              :: MK= KIND(1.0D0)
        INTEGER                         :: n_arg
        INTEGER                         :: num_dim
        INTEGER                         :: num_start
        INTEGER                         :: num_end
        INTEGER                         :: num_freq
        CHARACTER(256)                  :: cbuf
        INTEGER                         :: clen
        LOGICAL                         :: lExist
        LOGICAL                         :: lOpened
        CHARACTER(256)                  :: filename
        INTEGER                         :: fileunit
        INTEGER                         :: flength
        CHARACTER(256)                  :: unfoldname
        INTEGER                         :: unfoldunit
        REAL(MK),DIMENSION(:),POINTER   :: time
        REAL(MK),DIMENSION(3)           :: r
        REAL(MK),DIMENSION(3)           :: v 
        REAL(MK),DIMENSION(3)           :: f
        REAL(MK),DIMENSION(:),POINTER   :: theta
        REAL(MK),DIMENSION(:),POINTER   :: omega
        REAL(MK),DIMENSION(:),POINTER   :: torque
        
        INTEGER                         :: t
        REAL(MK)                        :: lx, rad
        INTEGER                         :: i,j
        INTEGER                         :: stat_info
        
        stat_info = 0
        num_start = 0
        num_end   = 1
        num_freq  = 1
        lx        = 5.0_MK
        rad       = 1.0_MK
        
        NULLIFY(time)
        NULLIFY(theta)
        NULLIFY(omega)
        NULLIFY(torque)
        
        n_arg = IARGC()
        
        IF (n_arg > 4 ) THEN
           
           CALL GETARG(1,cbuf)
           READ(cbuf,*), num_dim
           CALL GETARG(2,cbuf)
           READ(cbuf,*), num_start
           CALL GETARG(3,cbuf)
           READ(cbuf,*), num_end
           CALL GETARG(4,cbuf)
           READ(cbuf,*), num_freq
           CALL GETARG(5,cbuf)
           READ(cbuf,*), lx
           
        ELSE
           
           PRINT *, "command line argument not enough !"
           GOTO 9999
           
        END IF
        
        PRINT *, "num_dim : ", num_dim
        
        IF ( num_dim < 2 .OR. num_dim > 3 ) THEN
           PRINT *, "num_dim is not 2 or 3 ! "
           GOTO 9999
        END IF
        
        PRINT *, "num_start, num_end, num_freq :"
        PRINT *, num_start, num_end, num_freq
        
        IF ( num_start < 0 .OR. &
             num_end < 0 .OR. &
             num_end < num_start) THEN
           PRINT *, "num_start or num_end is wrong !"
           GOTO 9999
        END IF
        
        filename = TRIM("mcf_colloid01.dat")
        flength  = LEN_TRIM(filename)
        
        !------------------------------------
        ! Check if the particle file exists. 
        !------------------------------------
        
        lExist = .TRUE.
        INQUIRE(FILE=filename,EXIST=lExist)
        
        IF (.NOT.lExist) THEN
           WRITE(cbuf,'(2A)')'No such file: ', &
                filename(1:flength)
           PRINT *, cbuf
           GOTO 9999
        END IF
        
        stat_info = 0
        
        OPEN(fileunit,FILE=filename,&
             IOSTAT=stat_info,ACTION='READ')
        
        IF (stat_info /= 0) THEN
           WRITE(cbuf,'(2A)')'Failed to open file: ',&
                filename(1:flength)
           PRINT *, cbuf
           GOTO 9999
        END IF
        
        ALLOCATE(time(0:num_end))
        ALLOCATE(theta(0:num_end))
        ALLOCATE(omega(0:num_end))
        ALLOCATE(torque(0:num_end))

        !------------------------------------
        ! Scan the file and save data.
        !------------------------------------
        
        DO i = 0, num_end
           
           READ(fileunit,*,END=9999,ERR=200) &
                t,time(i), &
                r(1:num_dim), v(1:num_dim), &
                f(1:num_dim), theta(i), &
                omega(i), torque(i)
           
        END DO
        
        unfoldname =TRIM("mcf_colloid_r.dat")
        
        OPEN(UNIT=unfoldunit,FILE=unfoldname,&
             STATUS="REPLACE",FORM="FORMATTED",ACTION="WRITE", &
             IOSTAT=stat_info)
        
        IF (stat_info /= 0) THEN
           PRINT *, "Writting unfold file has problem !"
           GOTO 9999
        END IF
        

        WRITE(cbuf, '(A1,I2,A6)') '(', 2*num_dim + 1, 'E16.8)'
        clen = LEN_TRIM(cbuf)
        
        DO i = num_start, num_end, num_freq
           
           WRITE(unfoldunit,FMT=cbuf(1:clen),IOSTAT=stat_info) &
                time(i), theta(i), omega(i), torque(i)
           
        END DO
        
        PRINT *, "Unfold file of periodic image writtien to : ",&
             TRIM(unfoldname)
        
        GOTO 9999
        
200     CONTINUE
        
        PRINT *, "reading line ", i, " has error !"
        
        GOTO 9999
        
201     CONTINUE
        
        PRINT *, "reading the end already !"
        
9999    CONTINUE
        
        INQUIRE(UNIT=fileunit,OPENED=lOpened)
        IF (lopened) THEN
           CLOSE(fileunit)
        END IF
        
        INQUIRE(UNIT=unfoldunit,OPENED=lOpened)
        IF (lopened) THEN
           CLOSE(unfoldunit)
        END IF
        
        IF(ASSOCIATED(time)) THEN
           DEALLOCATE(time)
        END IF
        
        IF(ASSOCIATED(theta)) THEN
           DEALLOCATE(theta)
        END IF
        
        IF(ASSOCIATED(omega)) THEN
           DEALLOCATE(omega)
        END IF
        
        IF(ASSOCIATED(torque)) THEN
           DEALLOCATE(torque)
        END IF
        
      END PROGRAM unfold_r
      
      
