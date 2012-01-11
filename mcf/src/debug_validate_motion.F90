      SUBROUTINE debug_validate_motion_v(this,rank,caller,&
           step,vel,start,end,dt,cutoff,stat_info)
        !----------------------------------------------------
      	!  Subroutine   :        debug_validate_motion_v
      	!----------------------------------------------------
      	!
      	!  Purpose      : Subroutine of validate velocity and force
	!                 to ensure that movement of each time step
        !                 doesn't exceed cutoff length.
        !                 It handles 2 dimension data 
      	!
      	!  Input        : rank     (I) MPI rank of writing processor
      	!                 caller    (C) name of calling subroutine
	!		  dstep	    (I) step number
      	!                 dvel      (F) data pointer
	!                 dforce    (F) data point
    	!                 dt        (F) time step
	!
      	!  Input/output : 
      	!
      	!  Output       : stat_info      (I) error status
      	!
      	!  Routines     :
      	!
      	!  Remarks      :
      	!
      	!  References   :
      	!
      	!  Version	:
      	!----------------------------------------------------
      	! 
      	!----------------------------------------------------
        
      	!-------------------------------------------------------------
      	!	Arguments     
      	!-------------------------------------------------------------

        TYPE(Debug), INTENT(IN)                         :: this
        INTEGER, INTENT(IN)		                :: rank
        CHARACTER(LEN=*), INTENT(IN)		        :: caller
        INTEGER, INTENT(IN)	  			:: step
        REAL(MK), DIMENSION(:,:), POINTER 		:: vel
        INTEGER, INTENT(IN)                             :: start
        INTEGER, INTENT(IN)                             :: end
        REAL(MK), INTENT(IN)                            :: dt
	REAL(MK), INTENT(IN)                            :: cutoff
        INTEGER, INTENT(OUT)        	                :: stat_info

      	!-------------------------------------------------------------
      	!	Local variables 
      	!-------------------------------------------------------------

        INTEGER                         :: stat_info_sub
        INTEGER				:: i,j
        INTEGER, DIMENSION(2)		:: array_dim
        REAL(MK), DIMENSION(:),POINTER  :: dx
        INTEGER                         :: icaller
        
        REAL(MK)			:: time_routine_start
        INTEGER                         :: debug_threshold

      	!-------------------------------------------------------------
      	!	Initialize 
      	!-------------------------------------------------------------

        stat_info     = 0
	stat_info_sub = 0
        NULLIFY(dx)
        icaller = LEN_TRIM(caller)

        debug_threshold = 1

        IF( this%flag > debug_threshold) THEN
           CAll debug_substart(this,rank,&
                'debug_validate_motion', &
                time_routine_start,stat_info_sub)

        END IF
      

	IF ( start > end ) THEN

          CALL debug_print_msg(this,rank,&
               'debug_validate_motion_v_f',&
               "start index > end index happens !", stat_info_sub)
    
          stat_info = -1
          GOTO 9999	
          
       END IF
      

       array_dim(1) = SIZE(vel,1)
       array_dim(2) = SIZE(vel,2)
       ALLOCATE(dx(array_dim(1)))

       IF ( end - start + 1 > array_dim(2)) THEN

          CALL debug_print_msg(this,rank,&
               'debug_validate_motion_v', &
               'Velocity array exceeds bound !', stat_info_sub)
          
          stat_info = -1
          GOTO 9999
          
       END IF

    	!-------------------------------------------------------------
      	!	Validate the motiong
      	!------------------------------------------------------------- 

       DO j = start,end

          dx(1:array_dim(1)) = vel(1:array_dim(1),j) * dt

          DO i =1, array_dim(1)
             
             IF ( ABS(dx(i)) > cutoff) THEN
		 
                CALL debug_print_msg(this,rank, &
                     'debug_validate_motion_v',&
                     'One step movement exceeds cutoff',&
                     stat_info_sub) 
		 
                CALL debug_print_msg(this,rank, &
                     'debug_validate_motion_v',&
                     'Step',step,&
                     stat_info_sub) 

                CALL debug_print_msg(this,rank, &
                     'debug_validate_motion_v',&
                     'Particle index',j,&
                     stat_info_sub) 

                PRINT *, "vel:", vel(1:3,j)

                stat_info = -1
                GOTO 9999

             END IF

          END DO

       END DO

	!-------------------------------------------------------------
        !	Return 
      	!-------------------------------------------------------------


9999   CONTINUE

       IF( this%flag > debug_threshold) THEN

          CAll debug_substop(this,rank,&
               'debug_validate_motion_v', &
               time_routine_start,stat_info_sub)

       END IF
       
       
       RETURN
       
     END SUBROUTINE debug_validate_motion_v


	
      SUBROUTINE debug_validate_motion_v_f(this,rank,caller,&
           step,vel,force,start,end,dt,cutoff,stat_info)
        !-------------------------------------------------------------
      	!  Subroutine   :        debug_validate_motion_v_f
      	!-------------------------------------------------------------
      	!
      	!  Purpose      : Subroutine of validate velocity and force
	!                 to ensure that movement of each time step
        !                 doesn't exceed cutoff length.
        !                 It handles 2 dimension data 
      	!
      	!  Input        : rank     (I) MPI rank of writing processor
      	!                 caller    (C) name of calling subroutine
	!		  dstep	    (I) step number
      	!                 dvel      (F) data pointer
	!                 dforce    (F) data point
    	!                 dt        (F) time step
	!
      	!  Input/output : 
      	!
      	!  Output       : stat_info      (I) error status
      	!
      	!  Routines     :
      	!
      	!  Remarks      :
      	!
      	!  References   :
      	!
      	!  Version	:
      	!-------------------------------------------------------------
      	! 
      	!-------------------------------------------------------------
        
      	!-------------------------------------------------------------
      	!	Arguments     
      	!-------------------------------------------------------------
        TYPE(Debug), INTENT(IN)                         :: this
        INTEGER, INTENT(IN)		                :: rank
        CHARACTER(LEN=*), INTENT(IN)		        :: caller
        INTEGER, INTENT(IN)	  			:: step
        REAL(MK), DIMENSION(:,:), POINTER 		:: vel
        REAL(MK), DIMENSION(:,:), POINTER 		:: force
        INTEGER, INTENT(IN)                             :: start
        INTEGER, INTENT(IN)                             :: end
        REAL(MK), INTENT(IN)                            :: dt
	REAL(MK), INTENT(IN)                            :: cutoff
        INTEGER, INTENT(OUT)        	                :: stat_info

      	!-------------------------------------------------------------
      	!	Local variables 
      	!-------------------------------------------------------------

        INTEGER                         :: stat_info_sub
        INTEGER				:: i,j
        INTEGER, DIMENSION(2,2)		:: array_dim
        REAL(MK), DIMENSION(:),POINTER  :: dx
        INTEGER                         :: icaller

        REAL(MK)			:: time_routine_start
        INTEGER                         :: debug_threshold

      	!-------------------------------------------------------------
      	!	Initialize 
      	!-------------------------------------------------------------

        stat_info     = 0
	stat_info_sub = 0
        NULLIFY(dx)
        icaller = LEN_TRIM(caller)

        debug_threshold = 1

        IF( this%flag > debug_threshold) THEN

           CAll debug_substart(this,rank,&
                'debug_validate_motion', &
                time_routine_start,stat_info_sub)

        END IF
      

	IF ( start > end ) THEN

          CALL debug_print_msg(this,rank,&
               'debug_validate_motion_v_f',&
               "start index > end index happens !", stat_info_sub)
    
          stat_info = -1
          GOTO 9999	
          
       END IF
      

       array_dim(1,1) = SIZE(vel,1)
       array_dim(2,1) = SIZE(vel,2)
       ALLOCATE(dx(array_dim(1,1)))

       IF ( end - start + 1 > array_dim(2,1)) THEN

          CALL debug_print_msg(this,rank,&
               'debug_validate_motion_v_f', &
               'Velocity array exceeds bound !', stat_info_sub)
          
          stat_info = -1
          GOTO 9999
          
       END IF

       array_dim(1,2) = SIZE(force,1)
       array_dim(2,2) = SIZE(force,2)	

       IF ( end - start + 1 > array_dim(2,2)) THEN

          CALL debug_print_msg(this,rank,&
               'debug_validate_motion_v_f', &
               'Force array exceeds bound !', stat_info_sub)

          stat_info = -1
          GOTO 9999

       END IF

       IF (array_dim(1,1) /= array_dim(1,2) .OR. &
            array_dim(2,1) /= array_dim(2,2)) THEN

          CALL debug_print_msg(this,rank,&
               'debug_validate_motion_v_f', &
               'Velocity, Force arrays do not match !', stat_info_sub)

          stat_info = -1
          GOTO 9999

       END IF


	!-------------------------------------------------------------
      	!	Validate the motiong
      	!------------------------------------------------------------- 

       DO j = start,end

          dx(1:array_dim(1,1)) = vel(1:array_dim(1,1),j) * dt + &
               force(1:array_dim(1,1),j) * dt * dt

          DO i =1,array_dim(1,1)
             
             IF ( ABS(dx(i)) > cutoff) THEN
		 
                CALL debug_print_msg(this,rank, &
                     'debug_validate_motion_v_f',&
                     'One step movement exceeds cutoff',&
                     stat_info_sub) 
		 
                CALL debug_print_msg(this,rank, &
                     'debug_validate_motion_v_f',&
                     'Step',step,&
                     stat_info_sub) 

                CALL debug_print_msg(this,rank, &
                     'debug_validate_motion_v_f',&
                     'Particle index',j,&
                     stat_info_sub) 

                PRINT *, "vel:", vel(1:3,j)
                PRINT *, "force:", force(1:3,j)

                stat_info = -1
                GOTO 9999

             END IF

          END DO

       END DO

	!-------------------------------------------------------------
        !	Return 
      	!-------------------------------------------------------------


9999   CONTINUE


       IF( this%flag > debug_threshold) THEN
          
          CAll debug_substop(this,rank,&
               'debug_validate_motion_v_f', &
               time_routine_start,stat_info_sub)

       END IF
	   
	   
       RETURN
	
     END SUBROUTINE debug_validate_motion_v_f


