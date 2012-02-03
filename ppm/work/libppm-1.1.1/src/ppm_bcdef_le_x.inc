  !----------------------------------------------------------
  ! (Re)allocate memory for the Lees-Edwards ghosts.
  !----------------------------------------------------------
  
  iopt   = ppm_param_alloc_grow_preserve
  ldu(1) = ppm_dim
  ldu(2) = 2*nghostplus
  CALL ppm_alloc(xt,ldu,iopt,info) 
  IF (info.NE.0) THEN
     info = ppm_error_fatal
     CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
          'xt',__LINE__,info)
     GOTO 9999
  ENDIF
  
  ldu(1) = ldu(2)
  CALL ppm_alloc(ighost,ldu,iopt,info) 
  IF (info.NE.0) THEN
     info = ppm_error_fatal
     CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
          'ighost',__LINE__,info)
     GOTO 9999
  ENDIF
  
  !----------------------------------------------------------
  ! Lees-Edwards ghosts in the x-direction.
  !----------------------------------------------------------
  
  k = nghostplus
  
  !----------------------------------------------------------
  ! those at the west boundary.
  ! with '='.
  !----------------------------------------------------------
  
  xminf = min_phys(1) 
  xmini = min_phys(1) + ghostsize
  
  DO i = 1, nghostplus
     
     IF ( xt(1,i) >= xminf .AND. &
          xt(1,i) <  xmini ) THEN
        
        k         = k + 1
        ighost(k) = ighost(i)
        
        DO j = 1, ppm_dim
           
           IF ( j == 1 ) THEN
              
              xt(j,k) = xt(j,i) + len_phys(1)
              
           ELSE
              
              xt(j,k) = &
                   MODULO( xt(j,i) + shear_length(j,2), &
                   len_phys(j) )
              
           END IF
           
        END DO
        
     END IF
     
  END DO
  
  !----------------------------------------------------------
  ! east bc, if we are not using symmetry.
  ! with '='.
  !----------------------------------------------------------
  
  IF ( isymm == 0 ) THEN

     xmaxf = max_phys(1) 
     xmaxi = max_phys(1) - ghostsize
     
     
     DO i = 1, nghostplus
        
        IF ( xt(1,i) > xmaxi .AND. &
             xt(1,i) <= xmaxf ) THEN
           
           k         = k + 1
           ighost(k) = ighost(i)
           
           DO j = 1, ppm_dim
              
              IF ( j == 1 ) THEN
                 
                 xt(j,k) = xt(j,i) - len_phys(1)
                 
              ELSE
                 
                 xt(j,k) = &
                      MODULO(xt(j,i) + shear_length(j,1), &
                      len_phys(j) )
                 
              END IF
              
           END DO
           
        END IF
        
     END DO

  END IF ! isymm
  
  
  !----------------------------------------------------------
  ! update the ghost counter.
  !----------------------------------------------------------
  
  nghostplus = k
  
  
