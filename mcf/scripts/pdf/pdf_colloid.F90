      PROGRAM pdf_colloid
        !----------------------------------------------------
        ! Calculate Probability Distribution Fucntion of
        ! velocity for solvent particles and a colloid.
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
        INTEGER                         :: num_line_c
        REAL(MK)                        :: min_c
        REAL(MK)                        :: max_c
        INTEGER                         :: ndx_c
        REAL(MK)                        :: dx_c
        INTEGER                         :: num_c
        
        LOGICAL                         :: lExist
        LOGICAL                         :: lOpened
        
        CHARACTER(256)                  :: filein_c
        INTEGER                         :: unitin_c
        INTEGER                         :: lin_c
        
        CHARACTER(256)                  :: fileout_c
        INTEGER                         :: unitout_c
        INTEGER                         :: lout_c
        
        REAL(MK)                        :: time        
        REAL(MK),DIMENSION(3)           :: r
        REAL(MK),DIMENSION(3)           :: v
        
        INTEGER,DIMENSION(:),ALLOCATABLE:: pdf_c
        
        INTEGER                         :: i,j
        INTEGER                         :: ic
        
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
           READ(cbuf,*), num_line_c
           CALL GETARG(3,cbuf)
           READ(cbuf,*), max_c
           CALL GETARG(4,cbuf)
           READ(cbuf,*), ndx_c
           
        ELSE
           
           PRINT *, "command line argument not enough !"
           GOTO 9999
           
        END IF
        
        !----------------------------------------------------
        ! Define number of pieces and width for PDF.
        !----------------------------------------------------
        
        min_c = -max_c
        dx_c = (max_c-min_c) / ndx_c
        
        PRINT *, "num_dim    : ", num_dim
        PRINT *, "num_line_c : ", num_line_c
        PRINT *, "min_c      : ", min_c
        PRINT *, "max_c      : ", max_c
        PRINT *, "ndx_c      : ", ndx_c
        PRINT *, "dx_c       : ", dx_c
        
        !----------------------------------------------------
        ! Define file names
        !----------------------------------------------------
        
        filein_c = TRIM("mcf_colloid.dat")
        lin_c    = LEN_TRIM(filein_c)
        unitin_c  = 2
        
        fileout_c = TRIM("pdf_colloid.dat")
        lout_c    = LEN_TRIM(fileout_c)
        unitout_c = 4
   
        
        !----------------------------------------------------
        ! Check if the colloid file exists. 
        !----------------------------------------------------
        
        lExist = .TRUE.
        INQUIRE(FILE=filein_c,EXIST=lExist)
        
        IF (.NOT.lExist) THEN
           WRITE(cbuf,'(2A)')'No such file: ', &
                filein_c(1:lin_c)
           PRINT *, cbuf
           stat_info = -1           
           GOTO 9999
        END IF
        
        OPEN(unitin_c,FILE=filein_c,&
             IOSTAT=stat_info,ACTION='READ')
        
        IF (stat_info /= 0) THEN
           WRITE(cbuf,'(2A)')'Failed to open file: ',&
                filein_c(1:lin_c)
           PRINT *, cbuf
           stat_info = -1
           GOTO 9999
        END IF
        

        !----------------------------------------------------
        ! Scan the colloid file and save data.
        !----------------------------------------------------
        
        num_c = 0
        
        ALLOCATE(pdf_c(0:ndx_c))
        pdf_c(:) = 0
        
        DO i = 0, num_line_c
           
           READ(unitin_c,*,END=201,ERR=200) &
                time, r(1:num_dim), v(1:num_dim)
           
           DO j = 1, num_dim
              
              ic    = ANINT((v(j) - min_c) / dx_c)
              num_c = num_c + 1
              
              IF ( ic >= 0 .AND. ic <= ndx_c) THEN
                 
                 pdf_c(ic) = pdf_c(ic) + 1
                 
              ELSE
                 
                 PRINT *, "velocity of colloid goes out of given range !"
                 PRINT *, "v_c = ", v(j)
                 
              END IF
              
           END DO
           
        END DO
        
        !----------------------------------------------------
        ! Open output files.
        !----------------------------------------------------
        
        fileout_c =TRIM("pdf_colloid.dat")
        
        OPEN(UNIT=unitout_c,FILE=fileout_c,&
             STATUS="REPLACE",FORM="FORMATTED",ACTION="WRITE", &
             IOSTAT=stat_info)
        
        IF ( stat_info /= 0 ) THEN
           PRINT *, "Writting pdf_colloid file has problem !"
           stat_info = -1
           GOTO 9999
        END IF

        !----------------------------------------------------
        ! Writting pdf into files.
        !----------------------------------------------------
        
        
        DO i = 0, ndx_c
           
           WRITE(unitout_c,FMT='(2(E16.8))',IOSTAT=stat_info) &
                min_c + i * dx_c, 1.0_MK*pdf_c(i)/num_c/dx_c
           
        END DO
        
        PRINT *, "num_colloid : ", num_c
        PRINT *, "PDF of colloid velocity writtien to : ", &
             TRIM(fileout_c)
        
        
        GOTO 9999
        
        
200     CONTINUE
        
        PRINT *, "reading line ", i, " has error !"
        
        GOTO 9999
        
201     CONTINUE
        
        PRINT *, "reading the end already !"
        
9999    CONTINUE
        
        
        INQUIRE(UNIT=unitin_c,OPENED=lOpened)
        IF (lopened) THEN
           CLOSE(unitin_c)
        END IF
        
        INQUIRE(UNIT=unitout_c,OPENED=lOpened)
        IF (lopened) THEN
           CLOSE(unitout_c)
        END IF
        
        IF(ALLOCATED(pdf_c)) THEN
           DEALLOCATE(pdf_c)
        END IF
      
        
      END PROGRAM pdf_colloid
      
      
      
      
