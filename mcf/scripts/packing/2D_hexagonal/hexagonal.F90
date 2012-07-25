      PROGRAM hexagonal
        !----------------------------------------------------
        ! Generate positions of colloidal particles
        ! in 2D hexagonal pattern in a channel
        ! x-periodic. y-wall
        !----------------------------------------------------
        
        IMPLICIT NONE
        
        INTEGER                         :: n_arg
        CHARACTER(256)                  :: cbuf
        INTEGER, PARAMETER              :: MK= KIND(1.0D0)
        REAL(MK), DIMENSION(2)          :: box
        REAL(MK), DIMENSION(2)          :: dx
        REAL(MK), DIMENSION(2)          :: x_start
        REAL(MK), DIMENSION(2)          :: x
        
        REAL(MK)                        :: rad
        REAL(MK), DIMENSION(2)          :: gap
        REAL(MK)                        :: ratio1, ratio2
        INTEGER, DIMENSION(2)           :: num
        INTEGER                         :: flip
        INTEGER                         :: count
        INTEGER                         :: num_odd,num_even, num_t
        
        box(1)   = 32.0_MK
        box(2)   = 32.0_MK
        
        rad      = 1.0_MK
        gap(2)   = 1.0_MK
        
        n_arg = IARGC()

        IF (n_arg > 2 ) THEN
           
           CALL GETARG(1,cbuf)
           READ(cbuf,*), box(1)
           CALL GETARG(2,cbuf)
           READ(cbuf,*), box(2)
           CALL GETARG(3,cbuf)
           READ(cbuf,*), gap(2)
           
        END IF
        
     
        !num(2)   = INT(( box(2)/rad-ratio1) / ( 2 + ratio1))
        
        !----------------------------------------------------
        ! calculate number of partilces in y-direction
        !----------------------------------------------------
        
        num(2)   = INT(( box(2) - gap(2)) / ( 2.0_MK * rad + gap(2) ))
        
        !PRINT *, "num(2) : ", num(2)
        
        !----------------------------------------------------
        ! adjust gap(2)
        ! calculate dx(2) and x_start(2)
        !----------------------------------------------------

        gap(2) = (box(2) - 2.0_MK * num(2) * rad) / ( num(2) + 1)
        dx(2)  = 2.0_MK*rad+gap(2)
        
        dx(1)  = dx(2) * sqrt(3.0_MK)/ 2.0_MK
        gap(1) = dx(1)/2.0_MK
        
        num(1) = INT( box(1)/dx(1) )
        
        x(1)  = rad + dx(1)/2.0_MK 
        
        flip  = 1
        count = 0
        
        DO WHILE( x(1) <= box(1) - rad - gap(1) )
           
           flip = MOD(flip + 1, 2)
           
           x(2) = rad + gap(2) + dx(2) * flip / 2.0_MK
           
           IF ( flip==0 ) THEN
              num_even = 1
           ELSE
              num_odd = 1
           END IF
           
           DO WHILE ( x(2) < box(2) - rad )
              
              PRINT *, 'coll_shape = 1'
              PRINT *, 'coll_radius = ', rad,',', 0.0
              PRINT *, 'coll_x      = ', x(1),',', x(2)
              PRINT *, ' '
              count = count + 1
              x(2) = x(2) + dx(2)

              IF ( flip==0 ) THEN
                 num_even = num_even+1
              ELSE
                 num_odd = num_odd+1
              END IF

              IF ( num_odd == num_even -1 ) THEN
                 EXIT
              END IF
              
           END DO
           
           x(1) = x(1) + dx(1)
           
        END DO
        
        !PRINT *, box(1)
        PRINT *, ""
        PRINT *, "total number : ", count
        PRINT *, "volume fraction : ", &
             count * 3.141593_MK*rad**2/box(1)/box(2)
        
      END PROGRAM hexagonal
