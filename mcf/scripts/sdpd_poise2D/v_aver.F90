      PROGRAM v_aver
        
        IMPLICIT NONE
        
        INTEGER, PARAMETER              :: MK= KIND(1.0D0)
        
        INTEGER                         :: n_arg
        INTEGER                         :: num_line
        
        CHARACTER(256)                  :: cbuf
        LOGICAL                         :: lExist
        LOGICAL                         :: lOpened
        
        CHARACTER(256)                  :: filein_s
        INTEGER                         :: unitin_s
        INTEGER                         :: lin_s
           
        CHARACTER(256)                  :: fileout_s
        INTEGER                         :: unitout_s
        INTEGER                         :: lout_s
         
        REAL(MK)                        :: min_y
        REAL(MK)                        :: max_y
        INTEGER                         :: ndy
        REAL(MK)                        :: dy
        
        
        REAL(MK),DIMENSION(:,:),ALLOCATABLE :: v
        INTEGER,DIMENSION(:),ALLOCATABLE    :: num
        INTEGER                             :: index
        INTEGER                             :: i
        
        
        REAL(MK),DIMENSION(2)           :: r_s
        REAL(MK),DIMENSION(2)           :: v_s
        REAL(MK)                        :: rho
        REAL(MK)                        :: mass
        REAL(MK)                        :: pid
        REAL(MK)                        :: sid
        
        
        REAL(MK),PARAMETER              :: machine_zero=1.0e-10
        INTEGER                         :: stat_info
        
        !------------------------------------------
        ! Initialization
        !------------------------------------------
        
        stat_info = 0
        
        n_arg = IARGC()
        
        IF (n_arg > 1 ) THEN
           CALL GETARG(1,cbuf)
           READ(cbuf,*), num_line
           CALL GETARG(2,cbuf)
           READ(cbuf,*), ndy           
        ELSE
           PRINT *, "command line argument not enough !"
           GOTO 9999
        END IF
        
        PRINT *, "num_line, ndy : ",  num_line, ndy
        

        !--------------------------------
        ! Define number of pieces and
        ! width for PDF.
        !--------------------------------
        
        max_y = 1.0_MK
        min_y = 0.0_MK
        
        dy  = (max_y-min_y) / ndy
        
        ALLOCATE(v(2,ndy))
        ALLOCATE(num(ndy))
        
        v(:,:)   = 0.0_MK
        num(:) = 0
        
        !--------------------------------
        ! Define file names
        !--------------------------------
        
        filein_s = TRIM("particles.dat")
        lin_s    = LEN_TRIM(filein_s)
        unitin_s = 1
        
        fileout_s = TRIM("v.dat")
        lout_s    = LEN_TRIM(fileout_s)
        unitout_s = 3
        
        !------------------------------------
        ! Check if the particle file exists. 
        !------------------------------------
        
        lExist = .TRUE.
        INQUIRE(FILE=filein_s,EXIST=lExist)
        
        IF (.NOT.lExist) THEN
           WRITE(cbuf,'(2A)')'No such file: ', &
                filein_s(1:lin_s)
           PRINT *, cbuf
           stat_info = -1           
           GOTO 9999
        END IF
        
        OPEN(unitin_s,FILE=filein_s,&
             IOSTAT=stat_info,ACTION='READ')
        
        IF (stat_info /= 0) THEN
           WRITE(cbuf,'(2A)')'Failed to open file: ',&
                filein_s(1:lin_s)
           PRINT *, cbuf
           stat_info = -1
           GOTO 9999
        END IF
             
        !------------------------------------------
        ! Scan the particle file and save data.
        !------------------------------------------
        
        DO i = 1, num_line
           
           READ(unitin_s,*,END=201,ERR=200) &
                r_s(1:2), v_s(1:2), rho, mass, pid, sid
           
           IF (r_s(2) > max_y ) THEN
              
              r_s(2) = 2.0_MK * max_y - r_s(2)
              v_s(1) = -v_s(1)
              
           END IF
           
           index = INT(r_s(2)/dy) + 1
           
           IF ( index >0 .AND. index <= ndy   ) THEN
              
              v(1:2,index) = v(1:2,index) + v_s(1:2)
              num(index)   = num(index) + 1
              
           END IF
           
        END DO

        
        !------------------------------------------
        ! Open output file
        !------------------------------------------
        
        
        OPEN(UNIT=unitout_s,FILE=fileout_s,&
             STATUS="REPLACE",FORM="FORMATTED",ACTION="WRITE", &
             IOSTAT=stat_info)
        
        IF (stat_info /= 0) THEN
           PRINT *, "Writting v file has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !------------------------------------------
        ! Writting pdf into files
        !------------------------------------------
        
        DO i = 1, ndy
           
           r_s(2) = i * dy - 0.5_MK * dy
           
           WRITE(unitout_s,FMT='(3(E16.8))',IOSTAT=stat_info) &
                r_s(2), v(1:2,i)/num(index)
           
        END DO
        
        PRINT *, "Average velocity of solvent writtien to : ", TRIM(fileout_s)
        
        
        GOTO 9999
        
        
200     CONTINUE
        
        PRINT *, "reading line ", i, " has error !"
        
        GOTO 9999
        
201     CONTINUE
        
        PRINT *, "reading the end already !"
        
9999    CONTINUE
        
        INQUIRE(UNIT=unitin_s,OPENED=lOpened)
        IF (lopened) THEN
           CLOSE(unitin_s)
        END IF
        
        INQUIRE(UNIT=unitout_s,OPENED=lOpened)
        IF (lopened) THEN
           CLOSE(unitout_s)
        END IF
        
        IF(ALLOCATED(v)) THEN
           DEALLOCATE(v)
        END IF
        
        IF(ALLOCATED(num)) THEN
           DEALLOCATE(num)
        END IF
      
        
      END PROGRAM v_aver
      
      
      
      
