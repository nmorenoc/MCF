      SUBROUTINE debug_write_output_1d_f(this,rank,caller,prefix,&
           step,output,start,end,stat_info)
	!-------------------------------------------------------------
      	!  Subroutine :               debug_write_output
      	!-------------------------------------------------------------
      	!
      	!  Purpose    : Subroutine of writing data into file
      	!	   	for debugging purposes.
        !               It handles 1 dimension data 
      	!
      	!  Input      : rank     (I) MPI rank of writing processor
      	!               caller    (C) name of calling subroutine
	! 	        prefix    (A) prefix for the output name
	!		step	    (I) step number
      	!               output   (F) data 
	!		start    (I) start index
	!		end	    (I)	end index
      	!
      	! Input/output: 
      	!
      	!  Output     : stat_info      (I) error status
      	!
      	!  Routines   :
      	!
      	!  Remarks    :
      	!
      	!  References :
      	!
      	!  Version    :
      	!-------------------------------------------------------------
      	! 
      	!-------------------------------------------------------------

      	!-------------------------------------------------------------
      	!  Arguments
      	!-------------------------------------------------------------
        
        TYPE(Debug), INTENT(IN)                         :: this
        INTEGER,		INTENT(IN)		:: rank
        CHARACTER(LEN=*),	INTENT(IN)		:: caller
        CHARACTER(LEN=*),	INTENT(IN)		:: prefix
        INTEGER	  		  			:: step
        REAL(MK), DIMENSION(:),INTENT(IN)	        :: output
        INTEGER	  		  			:: start
        INTEGER	  		  			:: end 
        INTEGER                                         :: icaller
        INTEGER, INTENT(OUT)        	                :: stat_info
        
      	!-------------------------------------------------------------
      	! Local variables 
      	!-------------------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        CHARACTER(LEN=MAX_CHAR)		:: debug_filename
        INTEGER				:: debug_unit
        INTEGER				:: j
        INTEGER		                :: array_dim
        LOGICAL	 			:: lexist
        REAL(MK)			:: time_routine_start
        INTEGER                         :: debug_threshold
        
        !-------------------------------------------------------------
      	! Initialize 
      	!-------------------------------------------------------------
        
        stat_info = 0
        stat_info_sub = 0
        
        icaller = LEN_TRIM(caller)
        
        debug_threshold = 3
        
        IF ( this%flag > 2 .OR.&
             this%flag > debug_threshold) THEN
           CAll debug_substart(this,rank,'debug_write_output_1d', &
                time_routine_start,stat_info_sub)
        END IF
        
        IF ( start > end ) THEN
           CALL debug_print_msg(this,rank,'debug_write_output_1d_f',&
                "start index > end index happens !", stat_info_sub)
           stat_info = -1
           GOTO 9999
        END IF
        
        
        
        WRITE(debug_filename,'(A,I4.4,A,I8.8,A)') &
             'debug_output_'//prefix//'_rank', rank,'_', step,'.out'
        

        INQUIRE(FILE=debug_filename,EXIST=lexist)
        
        IF (lexist) THEN
           
           CALL debug_print_msg(this,rank,&
                'debug_write_output_1d',&
                "Existing debug file being overwritten !",&
                stat_info_sub)           
        END IF
        
        debug_unit = rank + 50
        OPEN(unit=debug_unit,file=debug_filename,&
             action="WRITE",status="REPLACE")
        
        array_dim = SIZE(output)
        IF ( this%flag > debug_threshold) THEN
           CALL debug_print_msg(this,rank,"debug_write_output_1d_f",&
                "array_dim",array_dim,stat_info_sub)
        END IF
        
        IF ( end - start > array_dim) THEN
           CALL debug_print_msg(this,rank,&
                'debug_write_output_1d_f', &
                'Excedes array bound !', stat_info_sub)
           PRINT *, "start, end, array_dim", start, end, array_dim
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( this%flag >debug_threshold) THEN
           CALL debug_print_msg(this,rank, &
                'debug_write_output_1d_f','array_dim',array_dim,stat_info_sub)
           CALL debug_print_msg(this,rank, &
                'debug_write_output_1d_f','start',start,stat_info_sub)
           CALL debug_print_msg(this,rank, &
                'debug_write_output_1d_f','end',end, stat_info_sub)
        END IF
        
        
	!-------------------------------------------------------------
      	! Do the writting debug output file
      	!-------------------------------------------------------------    
        
        
        DO j = start, end
	  
           WRITE(debug_unit,'(E16.8)'), output(j)	 

        END DO	       

	    
	 
        CLOSE(debug_unit)
	 
	!-------------------------------------------------------------
        ! Return 
      	!-------------------------------------------------------------
        

9999    CONTINUE

	
        IF ( this%flag > 2 .OR. this%flag > debug_threshold) THEN
    
           CAll debug_substop(this,rank,&
                'debug_write_output_1d_f', time_routine_start,stat_info_sub)
           
        END IF
        
        RETURN

      END SUBROUTINE debug_write_output_1d_f
      

      SUBROUTINE debug_write_output_1d_i(this,rank,caller,prefix,&
           step,output,start,end,stat_info)
        
      	!-------------------------------------------------------------
      	!	Arguments     
      	!-------------------------------------------------------------

        TYPE(Debug), INTENT(IN)                         :: this
      	INTEGER,		INTENT(IN)		:: rank
        CHARACTER(LEN=*),	INTENT(IN)		:: caller
        CHARACTER(LEN=*),	INTENT(IN)		:: prefix
        INTEGER	  		  			:: step
        INTEGER, DIMENSION(:), INTENT(IN)               :: output
        INTEGER	  		  			:: start
	INTEGER	  		  			:: end 
        INTEGER, INTENT(OUT)        	                :: stat_info

      	!-------------------------------------------------------------
      	!	Local variables 
      	!-------------------------------------------------------------

        INTEGER                         :: stat_info_sub
      	CHARACTER(LEN=MAX_CHAR)		:: debug_filename
        INTEGER				:: debug_unit
        INTEGER				:: j
        INTEGER		                :: array_dim
	LOGICAL	 			:: lexist
        REAL(MK)			:: time_routine_start
        INTEGER                         :: debug_threshold
        INTEGER                         :: icaller
        

      	!-------------------------------------------------------------
      	!	Initialize 
      	!-------------------------------------------------------------

	stat_info = 0
        stat_info_sub = 0
        icaller = LEN_TRIM(caller)
        
        debug_threshold = 3

        IF ( this%flag > 2 .OR.&
             this%flag > debug_threshold) THEN

	   CAll debug_substart(this,rank,'debug_write_output_1d_i', &
         time_routine_start,stat_info_sub)
    
        END IF

	IF ( start > end ) THEN
    
           CALL debug_print_msg(this,rank,'debug_write_output_1d_i',&
                "start index > end index happens !", stat_info_sub)

	   stat_info = -1
           GOTO 9999	
	   
	END IF
      
        WRITE(debug_filename,'(A,I4.4,A,I8.8,A)') &
              'debug_output_'//prefix//'_rank', rank,'_', step,'.out'


        INQUIRE(FILE=debug_filename,EXIST=lexist)

         IF (lexist) THEN

            CALL debug_print_msg(this,rank,&
                 'debug_write_output_1d_i',&
                 "Existing debug file being overwritten !",&
                 stat_info_sub)           
         END IF
	
         debug_unit = rank + 40
	
         OPEN(unit=debug_unit,file=debug_filename,&
              action="WRITE",status="REPLACE")

	array_dim = SIZE(output)

	
	IF ( this%flag > debug_threshold) THEN
		
	   CALL debug_print_msg(this,rank,"debug_write_output_1d_i",&
	   "array_dim",array_dim,stat_info_sub)
          
	END IF


	IF ( end - start > array_dim) THEN

	   CALL debug_print_msg(this,rank,&
	   'debug_write_output_1d_i', &
	   'Excedes array bound !', stat_info_sub)
            PRINT *, "start, end, array_dim", start, end, array_dim

	   stat_info = -1
	   GOTO 9999

	END IF


	IF ( this%flag >debug_threshold) THEN
	   
	   CALL debug_print_msg(this,rank, &
	   'debug_write_output_1d_i','array_dim',array_dim,stat_info_sub)	
	   CALL debug_print_msg(this,rank, &
	   'debug_write_output_1d_i','start',start,stat_info_sub)
	   CALL debug_print_msg(this,rank, &
	   'debug_write_output_1d_i','end',end, stat_info_sub)
	   
	END IF


	!-------------------------------------------------------------
      	!	Do the writting debug output file
      	!-------------------------------------------------------------    
 
        DO j = start, end
	  
           WRITE(debug_unit,'(I8)'), output(j)	 

        END DO	       

	    
	 
        CLOSE(debug_unit)
	 
	!-------------------------------------------------------------
        !	Return 
      	!-------------------------------------------------------------



9999    CONTINUE

	
        IF ( this%flag > 2 .OR. this%flag > debug_threshold) THEN
    
           CAll debug_substop(this,rank,&
                'debug_write_output_1d_f', time_routine_start,stat_info_sub)
           
        END IF
        
        RETURN

      END SUBROUTINE debug_write_output_1d_i
      

      SUBROUTINE debug_write_output_2d_f(this,rank,caller,prefix,&
           step,output,start,end,stat_info)
        !-------------------------------------------------------------
        ! Arguments
        !-------------------------------------------------------------
        
        TYPE(Debug), INTENT(IN)                         :: this
        INTEGER,		INTENT(IN)		:: rank
        CHARACTER(LEN=*),	INTENT(IN)		:: caller
        CHARACTER(LEN=*),	INTENT(IN)		:: prefix
        INTEGER	  		  			:: step
        REAL(MK), DIMENSION(:,:),INTENT(IN)             :: output
        INTEGER	  		  			:: start
        INTEGER	  		  			:: end
        INTEGER, INTENT(OUT)        	                :: stat_info
        
      	!-------------------------------------------------------------
      	! Local variables 
      	!-------------------------------------------------------------
        
        INTEGER                         :: stat_info_sub
      	CHARACTER(LEN=MAX_CHAR)		:: debug_filename
        INTEGER				:: debug_unit
        INTEGER				:: j
        INTEGER, DIMENSION(2)		:: array_dim
        LOGICAL	 			:: lexist
        REAL(MK)			:: time_routine_start
        INTEGER                         :: debug_threshold
        INTEGER                         :: icaller
        
      	!-------------------------------------------------------------
      	! Initialize 
      	!-------------------------------------------------------------
        
        stat_info = 0
        stat_info_sub = 0
        
        icaller = LEN_TRIM(caller)
        
        debug_threshold = 3

        IF ( this%flag > 2 .OR.&
             this%flag >debug_threshold) THEN
           CAll debug_substart(this,rank,'debug_write_output_2d_f', &
                time_routine_start,stat_info_sub)
        END IF
        
        IF ( start > end ) THEN
           
           CALL debug_print_msg(this,rank,'debug_write_output_2d_f',&
                "start index > end index happens !", stat_info_sub)
           stat_info = -1
           GOTO 9999	
           
        END IF
        
        WRITE (debug_filename,'(A,I4.4,A,I8.8,A)') &
             'debug_output_'//prefix//'_rank', rank,'_', step,'.out'
        
        
        INQUIRE(FILE=debug_filename,EXIST=lexist)
        
        IF (lexist) THEN
           
           CALL debug_print_msg(this,rank,&
                'debug_write_output_2d_f',&
                "Existing debug file being overwritten !",&
                stat_info_sub)           
        END IF
	
        debug_unit = rank + 40
	
        OPEN(unit=debug_unit,file=debug_filename,&
             action="WRITE",status="REPLACE")
        
        array_dim(1) = SIZE(output,1)
        array_dim(2) = SIZE(output,2)
        
        IF ( this%flag > debug_threshold ) THEN
           CALL debug_print_msg(this,rank,"debug_write_output_2d_f",&
                "array_dim(1)",array_dim(1),stat_info_sub)
           CALL debug_print_msg(this,rank,"debug_write_output_2d_f",&
                "array_dim(2)",array_dim(2),stat_info_sub)
           
        END IF

        IF ( end - start > array_dim(2)) THEN
           CALL debug_print_msg(this,rank,&
                'debug_write_output_2d_f', &
                'Excedes array bound !', stat_info_sub)
           PRINT *, "start, end, array_dim(2)", start, end, array_dim(2)
           stat_info = -1 
           GOTO 9999
        END IF
        
        IF ( this%flag >debug_threshold) THEN
           
           CALL debug_print_msg(this,rank, &
                'debug_write_output_2d_f','array_dim(1)',&
                array_dim(1),stat_info_sub)
           CALL debug_print_msg(this,rank, &
                'debug_write_output_2d_f','array_dim(2)',&
                array_dim(2),stat_info_sub)
           CALL debug_print_msg(this,rank, &
                'debug_write_output_2d_f','start',start,stat_info_sub)
           CALL debug_print_msg(this,rank, &
                'debug_write_output_2d_f','end',end, stat_info_sub)
        END IF
        
        !----------------------------------------------------
        ! do the writting of debug output file
      	!----------------------------------------------------

        IF ( array_dim(1) == 1) THEN
           
           DO j = start, end
              WRITE(debug_unit,'(E16.8)'), output(1,j)	 
           END DO
           
        ELSE IF ( array_dim(1) == 2) THEN
           
           DO j = start, end
              WRITE(debug_unit,'(2E16.8)'), output(1:2,j)
           END DO
           
        ELSE IF (  array_dim(1) == 3) THEN
           
           DO j = start, end
              WRITE(debug_unit,'(3E16.8)'), output(1:3,j)
           END DO
           
        ELSE IF (  array_dim(1) == 4) THEN
           
           DO j = start, end
              WRITE(debug_unit,'(4E16.8)'), output(1:4,j)
           END DO
           
        ELSE IF (  array_dim(1) == 5) THEN
           
           DO j = start, end
              WRITE(debug_unit,'(5E16.8)'), output(1:5,j)
           END DO

        ELSE IF (  array_dim(1) == 6) THEN
           
           DO j = start, end
              WRITE(debug_unit,'(6E16.8)'), output(1:6,j)
           END DO
           
        ELSE IF (  array_dim(1) == 7) THEN
           
           DO j = start, end
              WRITE(debug_unit,'(7E16.8)'), output(1:7,j)
           END DO
           
        ELSE IF (  array_dim(1) == 8) THEN
           
           DO j = start, end
              WRITE(debug_unit,'(8E16.8)'), output(1:8,j)
           END DO

        ELSE IF (  array_dim(1) == 9) THEN
           
           DO j = start,end
              WRITE(debug_unit,'(9E16.8)'), output(1:9,j)
           END DO

        ELSE

           CALL debug_print_msg(this,rank,&
                "debug_write_output_2d_f",&
                "Number of rows can't be handled!",stat_info_sub)
           stat_info = -1
           GOTO 9999
           
        END IF
        
        
        CLOSE(debug_unit)
        
	!----------------------------------------------------
        ! Return 
      	!----------------------------------------------------

        
        
9999    CONTINUE
        
        IF ( this%flag > 2 .OR. this%flag > debug_threshold) THEN
           
           CAll debug_substop(this,rank,&
                'debug_write_output_2d_f', time_routine_start,stat_info)
           
        END IF
        
        RETURN
        
      END SUBROUTINE debug_write_output_2d_f
      

      SUBROUTINE debug_write_output_2d_i(this,rank,caller,prefix,&
           step,output,start,end,stat_info)

      	!-------------------------------------------------------------
      	!	Arguments     
      	!-------------------------------------------------------------

        TYPE(Debug), INTENT(IN)                         :: this
      	INTEGER,		INTENT(IN)		:: rank
        CHARACTER(LEN=*),	INTENT(IN)		:: caller
        CHARACTER(LEN=*),	INTENT(IN)		:: prefix
        INTEGER	  		  			:: step
        INTEGER, DIMENSION(:,:), INTENT(IN) 		:: output
        INTEGER	  		  			:: start
	INTEGER	  		  			:: end 
        INTEGER, INTENT(OUT)        	                :: stat_info

      	!-------------------------------------------------------------
      	!	Local variables 
      	!-------------------------------------------------------------

        INTEGER                         :: stat_info_sub
      	CHARACTER(LEN=MAX_CHAR)		:: debug_filename
        INTEGER				:: debug_unit
        INTEGER				:: j
        INTEGER, DIMENSION(2)		:: array_dim
	LOGICAL	 			:: lexist
        REAL(MK)			:: time_routine_start
        INTEGER                         :: debug_threshold
        INTEGER                         :: icaller
        
      	!-------------------------------------------------------------
      	!	Initialize 
      	!-------------------------------------------------------------

	stat_info = 0
        stat_info_sub = 0
        icaller = LEN_TRIM(caller)
        
        debug_threshold = 3
        
        IF ( this%flag > 2 .OR.&
             this%flag >debug_threshold) THEN

	   CAll debug_substart(this,rank,'debug_write_output_2d_i', &
         time_routine_start,stat_info_sub)
    
        END IF

	IF ( start > end ) THEN
    
           CALL debug_print_msg(this,rank,'debug_write_output_2d_i',&
                "start index > end index happens !", stat_info_sub)

	   stat_info = -1
           GOTO 9999	
	   
	END IF
      
         WRITE(debug_filename,'(A,I4.4,A,I8.8,A)') &
              'debug_output_'//prefix//'_rank', rank,'_', step,'.out'


         INQUIRE(FILE=debug_filename,EXIST=lexist)

         IF (lexist) THEN

            CALL debug_print_msg(this,rank,&
                 'debug_write_output_2d_i',&
                 "Existing debug file being overwritten !",&
                 stat_info_sub)           
         END IF
	
	debug_unit = rank + 40
      	OPEN(unit=debug_unit,file=debug_filename,&
	action="WRITE",status="REPLACE")

	array_dim(1) = SIZE(output,1)
	array_dim(2) = SIZE(output,2)	

	
	IF ( this%flag > debug_threshold) THEN
		
	   CALL debug_print_msg(this,rank,"debug_write_output_2d_i",&
	   "array_dim(1)",array_dim(1),stat_info_sub)
           CALL debug_print_msg(this,rank,"debug_write_output_2d_i",&
	   "array_dim(2)",array_dim(2),stat_info_sub)
	   
	END IF


	IF ( end - start > array_dim(2)) THEN

	   CALL debug_print_msg(this,rank,&
	   'debug_write_output_2d_i', &
	   'Excedes array bound !', stat_info_sub)
           PRINT *, "start, end, array_dim(2)", start, end, array_dim(2)

	   stat_info = -1
	   GOTO 9999

	END IF


	IF ( this%flag >debug_threshold) THEN
	   
	   CALL debug_print_msg(this,rank, &
	   'debug_write_output_2d_i','array_dim(1)',array_dim(1),stat_info_sub)
	   CALL debug_print_msg(this,rank, &
	   'debug_write_output_2d_i','array_dim(2)',array_dim(2),stat_info_sub)
	   CALL debug_print_msg(this,rank, &
	   'debug_write_output_2d_i','start',start,stat_info_sub)
	   CALL debug_print_msg(this,rank, &
	   'debug_write_output_2d_i','end',end, stat_info_sub)
	   
	END IF


	!-------------------------------------------------------------
      	!	Do the writting debug output file
      	!-------------------------------------------------------------    

	    
	IF ( array_dim(1) == 1) THEN
	       
	   DO j = start, end
		  
	      WRITE(debug_unit,'(I8)'), output(1,j)	 

	   END DO
	       
	ELSE IF ( array_dim(1) == 2) THEN
	     
	   DO j = start, end
  
	      WRITE(debug_unit,'(2I8)'), output(1:2,j)

	   END DO
	   
	ELSE IF (  array_dim(1) == 3) THEN

	   DO j = start, end
	      
	      WRITE(debug_unit,'(3I8)'), output(1:3,j)	 
	      
	   END DO

	ELSE IF (  array_dim(1) == 4) THEN
	   
	   DO j = start, end

	      WRITE(debug_unit,'(4I8)'), output(1:4,j)
	      
	   END DO

	ELSE IF (  array_dim(1) == 5) THEN

	   DO j = start, end

	      WRITE(debug_unit,'(5I8)'), output(1:5,j)
		      
	   END DO

	ELSE IF (  array_dim(1) == 6) THEN
	       
	   DO j = start, end

	      WRITE(debug_unit,'(6I8)'), output(1:6,j) 
		
	   END DO
	       
	ELSE IF (  array_dim(1) == 7) THEN

	   DO j = start, end
	       
	      WRITE(debug_unit,'(7I8)'), output(1:7,j)
 
	   END DO
	   
	ELSE IF (  array_dim(1) == 8) THEN

	   DO j = start, end

	      WRITE(debug_unit,'(8I8)'), output(1:8,j)

	   END DO
	
	ELSE IF (  array_dim(1) == 9) THEN
	       
	   DO j = start,end

	      WRITE(debug_unit,'(9I8)'), output(1:9,j)

	   END DO

       ELSE 
	       
	   CALL debug_print_msg(this,rank,&
	   "debug_write_output_2d_i",&
	   "Number of rows can't be handled!",stat_info_sub)
	       
	   stat_info = -1
	   GOTO 9999
	       
        END IF
	    
	 
        CLOSE(debug_unit)
	 
	!-------------------------------------------------------------
        !	Return 
      	!-------------------------------------------------------------



9999    CONTINUE

	
        IF ( this%flag > 2 .OR. this%flag > debug_threshold) THEN
    
           CAll debug_substop(this,rank,&
                'debug_write_output_2d_i', time_routine_start,stat_info)
           
        END IF
        
        RETURN

      END SUBROUTINE debug_write_output_2d_i
      

