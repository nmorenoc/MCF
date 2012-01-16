      SUBROUTINE debug_print_msg_a(this,rank,caller,cbuf,stat_info)
        !----------------------------------------------------
      	! Subroutine  : debug_print_msg_a
      	!----------------------------------------------------
	!
      	! Purpose     : Subroutine of printing string message 
      	!	 	for debugging purposes.
      	!
      	! Input       : rank   (I) MPI rank of writing processor
      	!               caller (C) name of calling subroutine
      	!               cbuf   (C) message
      	!
      	!Input/output : 
      	!
      	! Output      : stat_info      (I) error status
      	!
      	! Routines    :
      	!
      	! Remarks     :
      	!
      	! References  :
      	!
      	! Revisions   : V0.1 01.03.2009, original version.
	!
      	!----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
      
        
      	!----------------------------------------------------
      	! Arguments     
      	!----------------------------------------------------

        TYPE(Debug), INTENT(IN)         :: this
        INTEGER			        :: rank
      	CHARACTER(LEN=*) 	        :: caller,cbuf
        INTEGER          	        :: stat_info
       
      	!----------------------------------------------------
      	! Local variables 
      	!----------------------------------------------------

        CHARACTER(LEN=MAX_CHAR)	        :: cformat, ctitle
        INTEGER            	        :: ios
        INTEGER            	        :: icaller,ibuf
        
      	!----------------------------------------------------
      	! Initialize 
      	!----------------------------------------------------
        
        stat_info = this%flag
        stat_info = 0
        
      	!----------------------------------------------------
      	! Get length of the caller's name and messages
      	!----------------------------------------------------
        
        icaller = LEN_TRIM(caller)
        ibuf    = LEN_TRIM(cbuf)
        
      	!----------------------------------------------------
      	! Define the print format
        !----------------------------------------------------
        
        IF (rank < -1) THEN
           
           cformat = '(3A)'
           
        ELSE IF (rank < 0) THEN
           
           cformat = '(A,I2,3A)'
           
        ELSE IF (rank.LT.10) THEN
           
           cformat = '(A,I1,3A)' 
           
        ELSE IF (rank.LT.100) THEN
           
           cformat = '(A,I2,3A)' 
           
        ELSE IF (rank.LT.1000) THEN
           
           cformat = '(A,I3,3A)' 
        ELSE
           
           cformat = '(A,I4,3A)' 
           
        END IF
        
        !----------------------------------------------------
        ! Do the printing
        !----------------------------------------------------
        
        IF (rank < -1) THEN
           
           WRITE(ctitle,cformat,IOSTAT=ios) &
                &      '<',                 &
                &      caller(1:icaller),   &
                &      '> : '
           
        ELSE
           
           WRITE(ctitle,cformat,IOSTAT=ios) &
                &      '[',rank,'] <',      &
                &      caller(1:icaller),   &
                &      '> : '
           
        END IF
        
        WRITE(*, '(A45A)',IOSTAT=ios), ctitle, cbuf(1:ibuf)
        
        !----------------------------------------------------
        ! Return 
      	!----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE debug_print_msg_a
      
      
      
      SUBROUTINE debug_print_msg_i(this,rank,caller,inum,stat_info)
        !----------------------------------------------------
      	! Subroutine  :  debug_print_msg_i
      	!----------------------------------------------------
	!
      	! Purpose     : Subroutine of printing integer message 
      	!	 	for debugging purposes.
      	!
      	! Input       : rank     (I) MPI rank of writing processor
      	!               caller    (C) name of calling subroutine
      	!               cbuf      (C) message
      	!
      	! Input/output: 
      	!
      	! Output      : stat_info      (I) error status
      	!
      	! Routines    :
      	!
      	! Remarks     :
      	!
      	! References  :
      	!
      	! Revisions   : V0.1 01.03.2009, original version.
	!
      	!----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
      

      	!----------------------------------------------------
      	! Arguments     
      	!----------------------------------------------------
        
        TYPE(Debug), INTENT(IN)         :: this
        INTEGER			        :: rank
        CHARACTER(LEN=*) 	        :: caller
        INTEGER			        :: inum
        INTEGER          	        :: stat_info
        
        !----------------------------------------------------
        ! Local variables 
      	!----------------------------------------------------
        
        CHARACTER(LEN=MAX_CHAR)	:: ibuf
        
        !----------------------------------------------------
        ! Initialize 
      	!----------------------------------------------------
        
        stat_info    = 0
        
        WRITE(ibuf,'(I10)'), inum
        
        CALL debug_print_msg_a(this,rank,caller,ibuf,stat_info)
        
      	!----------------------------------------------------
      	! Return 
      	!----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE debug_print_msg_i
      
      
      SUBROUTINE debug_print_msg_f(this,rank,caller,fnum,stat_info)
        !----------------------------------------------------
        ! Subroutine   :         debug_print_msg_f
      	!----------------------------------------------------
	!
      	! Purpose      : Subroutine of printing real number message 
      	!   	      	  for debugging purposes.
      	!
      	! Input        : rank     (I) MPI rank of writing processor
      	!                 caller    (C) name of calling subroutine
      	!                 fnum      (C) message
      	!
      	! Input/output : 
      	!
      	! Output       : stat_info      (I) error status
      	!
      	! Routines     :
      	!
      	! Remarks      :
      	!
      	! References   :
        !
     	! Revisions   : V0.1 01.03.2009, original version.
	!
      	!----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
         	
        
      	!----------------------------------------------------
      	! Arguments     
      	!----------------------------------------------------
        
        TYPE(Debug), INTENT(IN)         :: this
        INTEGER			        :: rank
        CHARACTER(LEN=*) 	        :: caller
        REAL(MK)		        :: fnum
        INTEGER          	        :: stat_info
        
      	!----------------------------------------------------
      	! Local variables 
      	!----------------------------------------------------
        
        CHARACTER(LEN=MAX_CHAR)	:: fbuf
        
        !----------------------------------------------------
        ! Initialize 
      	!----------------------------------------------------
        
        stat_info    = 0
        
        WRITE(fbuf,'(F11.8)'), fnum
        
        CALL debug_print_msg_a(this,rank,caller,fbuf,stat_info)
        
        !----------------------------------------------------
      	! Return 
      	!----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE debug_print_msg_f
      
      
      SUBROUTINE debug_print_msg_aa(this,rank,caller,abuf,bbuf,stat_info)
        !-------------------------------------------------------------
      	! Subroutine   :        debug_print_msg_aa
      	!-------------------------------------------------------------
	!
      	! Purpose      : Subroutine of printing string message 
      	!	 	      	for debugging purposes.
      	!
      	! Input        : rank     (I) MPI rank of writing processor
      	!                 caller    (C) name of calling subroutine
	!		  abuf	    (A) message
      	!                 bbuf      (A) string message
      	!
      	! Input/output : 
      	!
      	! Output       : stat_info      (I) error status
      	!
      	! Routines     :
      	!
      	! Remarks      :
      	!
      	! References   :
      	!
        ! Revisions   : V0.1 01.03.2009, original version.
	!
      	!----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
      
        
      	!----------------------------------------------------
      	! Arguments     
      	!----------------------------------------------------
        
        TYPE(Debug), INTENT(IN)         :: this
        INTEGER			        :: rank
        CHARACTER(LEN=*) 	        :: caller,abuf,bbuf
        INTEGER          	        :: stat_info
        
        !----------------------------------------------------
      	! Local variables 
      	!----------------------------------------------------
        
        CHARACTER(LEN=MAX_CHAR)	:: cbuf
        
        !----------------------------------------------------
        ! Initialize 
      	!----------------------------------------------------
        
        stat_info  = 0
       
        WRITE(cbuf,'(A,A)'), abuf//' : ',bbuf
        
        CALL debug_print_msg_a(this,rank,caller,cbuf,stat_info)
        
        !----------------------------------------------------
      	! Return 
      	!----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE debug_print_msg_aa
      
      
      SUBROUTINE debug_print_msg_ai(this,rank,caller,abuf,inum,stat_info)
        !----------------------------------------------------
        ! Subroutine   : debug_print_msg_ai
      	!----------------------------------------------------
	!
      	! Purpose      : Subroutine of printing integer message 
      	!	 	      	for debugging purposes.
      	!
      	! Input        : rank     (I) MPI rank of writing processor
      	!                 caller    (C) name of calling subroutine
	!		  abuf	    (A) message
      	!                 inum      (C) integer message
      	!
      	! Input/output : 
      	!
      	! Output       : stat_info      (I) error status
      	!
      	! Routines     :
      	!
      	! Remarks      :
      	!
      	! References   :
       	! 
        ! Revisions   : V0.1 01.03.2009, original version.
	!
      	!----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
     

      	!----------------------------------------------------
      	! Arguments     
      	!----------------------------------------------------
        
        TYPE(Debug), INTENT(IN)         :: this
        INTEGER			        :: rank
        CHARACTER(LEN=*) 	        :: caller,abuf
        INTEGER			        :: inum
        INTEGER          	        :: stat_info
        
        !----------------------------------------------------
      	! Local variables 
      	!----------------------------------------------------
        
        CHARACTER(LEN=MAX_CHAR)	:: cbuf
        
        !----------------------------------------------------
        ! Initialize 
      	!----------------------------------------------------
        
        stat_info    = 0
        
        WRITE(cbuf,'(A,I10)'), abuf//' : ',inum
        
        CALL debug_print_msg_a(this,rank,caller,cbuf,stat_info)
        
	
        !----------------------------------------------------
      	! Return 
      	!----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE debug_print_msg_ai
      
      
      SUBROUTINE debug_print_msg_af(this,rank,caller,abuf,fnum,stat_info)
        !----------------------------------------------------
      	! Subroutine  :  debug_print_msg_af
      	!----------------------------------------------------
	!
      	! Purpose     : Subroutine of printing integer message 
      	! 	      	for debugging purposes.
      	!
      	! Input       : rank    (I) MPI rank of writing processor
      	!               caller  (C) name of calling subroutine
	!               abuf	(A) message
      	!               fnum    (C) real number message
      	!
      	! Input/output: 
      	!
      	! Output      : stat_info      (I) error status
      	!
      	! Routines    :
      	!
      	! Remarks     :
      	!
      	! References  :
      	!
        ! Revisions   : V0.1 01.03.2009, original version.
	!
      	!----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
     

      	!----------------------------------------------------
      	! Arguments     
      	!----------------------------------------------------
        
        TYPE(Debug), INTENT(IN)         :: this
        INTEGER			        :: rank
        CHARACTER(LEN=*) 	        :: caller,abuf	
        REAL(MK)		        :: fnum
        INTEGER          	        :: stat_info
       
        !----------------------------------------------------
      	! Local variables 
      	!----------------------------------------------------
        
        CHARACTER(LEN=MAX_CHAR)	:: cbuf
        
        !----------------------------------------------------
      	! Initialize 
      	!----------------------------------------------------
        
        stat_info    = 0
        WRITE(cbuf,'(A,F11.8)'), abuf//' : ',fnum
        
        CALL debug_print_msg_a(this,rank,caller,cbuf,stat_info)
        
      	!----------------------------------------------------
      	! Return 
      	!----------------------------------------------------
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE debug_print_msg_af
