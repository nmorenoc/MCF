      PROGRAM pdf_solvent
        !----------------------------------------------------
        ! Calculate Probability Distribution Fucntion of
        ! velocity for solvent particles
        !
        ! Input 
        !
        ! num_dim    : number of dimension
        ! num_line_s : number of lines in particles file.
        ! num_line_c : number of lines in colloid file.
        !
        ! Output   :
        !             v, pdf
        !----------------------------------------------------
        
        IMPLICIT NONE
        
        INTEGER, PARAMETER              :: MK= KIND(1.0D0)
        
        INTEGER                         :: n_arg
        CHARACTER(256)                  :: cbuf
        INTEGER                         :: num_dim
        INTEGER                         :: num_line_s
        REAL(MK)                        :: min_s
        REAL(MK)                        :: max_s
        INTEGER                         :: ndx_s
        REAL(MK)                        :: dx_s
        INTEGER                         :: num_s

        LOGICAL                         :: lExist
        LOGICAL                         :: lOpened
        
        CHARACTER(256)                  :: filein_s
        INTEGER                         :: unitin_s
        INTEGER                         :: lin_s
        
        CHARACTER(256)                  :: fileout_s
        INTEGER                         :: unitout_s
        INTEGER                         :: lout_s
        
        
        
        REAL(MK),DIMENSION(3)           :: r
        REAL(MK),DIMENSION(3)           :: v
        REAL(MK)                        :: rho
        REAL(MK)                        :: mass
        REAL(MK)                        :: pid
        REAL(MK)                        :: sid
        
        INTEGER,DIMENSION(:),ALLOCATABLE:: pdf_s
        
        INTEGER                         :: i,j
        INTEGER                         :: is
        
        REAL(MK),PARAMETER              :: machine_zero=1.0e-6
        INTEGER                         :: stat_info
        
        !----------------------------------------------------
        ! Initialization
        !----------------------------------------------------
        
        stat_info = 0
        
        n_arg = IARGC()
        
        IF ( n_arg > 3 ) THEN
           
           CALL GETARG(1,cbuf)
           READ(cbuf,*), num_dim
           CALL GETARG(2,cbuf)
           READ(cbuf,*), num_line_s
           CALL GETARG(3,cbuf)
           READ(cbuf,*), max_s
           CALL GETARG(4,cbuf)
           READ(cbuf,*), ndx_s
           
        ELSE
           
           PRINT *, "command line argument not enough !"
           GOTO 9999
           
        END IF
        
        !----------------------------------------------------
        ! Define number of pieces and width for PDF.
        !----------------------------------------------------
        
        min_s = -max_s
        dx_s = (max_s-min_s) / ndx_s
        
        PRINT *, "num_dim     : ", num_dim
        PRINT *, "num_line_s  : ", num_line_s
        PRINT *, "min_s       : ", min_s
        PRINT *, "max_s       : ", max_s
        PRINT *, "ndx_s       : ", ndx_s
        PRINT *, "dx_s        : ", dx_s
        
        
        !----------------------------------------------------
        ! Define file names
        !----------------------------------------------------
        
        filein_s = TRIM("particles.dat")
        lin_s    = LEN_TRIM(filein_s)
        unitin_s = 1
        
        fileout_s = TRIM("pdf_solvent.dat")
        lout_s    = LEN_TRIM(fileout_s)
        unitout_s = 3
        
        !----------------------------------------------------
        ! Check if the particle file exists. 
        !----------------------------------------------------
        
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
        
        
        !----------------------------------------------------
        ! Scan the particle file and save data.
        !----------------------------------------------------
        
        num_s = 0
        
        ALLOCATE(pdf_s(0:ndx_s))
        pdf_s(:) = 0
        
        DO i = 0, num_line_s
           
           READ(unitin_s,*,END=201,ERR=200) &
                r(1:num_dim), v(1:num_dim), rho, mass, pid,sid
           
           IF ( ABS(sid) < machine_zero ) THEN
              
              DO j = 1, num_dim
                 
                 is = ANINT((v(j) - min_s) / dx_s)
                 num_s = num_s + 1
                 
                 IF ( is >= 0 .AND. is <= ndx_s) THEN
                    
                    pdf_s(is) = pdf_s(is) + 1
                    
                 ELSE
                    
                    PRINT *, "velocity of particle goes out of given range !"
                    PRINT *, "v_s = ", v(j)
                    
                 END IF
                 
              END DO
              
           END IF
           
        END DO

        !----------------------------------------------------
        ! Open output files.
        !----------------------------------------------------
        
        fileout_s =TRIM("pdf_solvent.dat")
        
        OPEN(UNIT=unitout_s,FILE=fileout_s,&
             STATUS="REPLACE",FORM="FORMATTED",ACTION="WRITE", &
             IOSTAT=stat_info)
        
        IF ( stat_info /= 0 ) THEN
           PRINT *, "Writting pdf_solvent file has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Writting pdf into files.
        !----------------------------------------------------
        
        DO i = 0, ndx_s
           
           WRITE(unitout_s,FMT='(2(E16.8))',IOSTAT=stat_info) &
                min_s + i * dx_s, 1.0_MK*pdf_s(i)/num_s/dx_s
           
        END DO
        
        PRINT *, "num_solvent : ", num_s
        PRINT *, "PDF of solvent velocity writtien to : ", &
             TRIM(fileout_s)
        
        
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
        
        IF(ALLOCATED(pdf_s)) THEN
           DEALLOCATE(pdf_s)
        END IF
        
        
      END PROGRAM pdf_solvent
      
      
      
      
