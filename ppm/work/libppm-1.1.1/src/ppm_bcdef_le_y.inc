  !----------------------------------------------------------
  ! (Re)allocate memory for the Leeds-Edwards ghosts.
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
  ! Lees-Edwards ghosts in the y-direction
  !----------------------------------------------------------
  
  k  = nghostplus
  
  !----------------------------------------------------------
  ! those at the south boundary 
  ! wiht '=';
  !----------------------------------------------------------
  
  yminf = min_phys(2) 
  ymini = min_phys(2) + ghostsize
  
  DO i = 1, nghostplus
     
     IF ( xt(2,i) >= yminf .AND. &
          xt(2,i) <  ymini ) THEN
        
        k         = k + 1
        ighost(k) = ighost(i)
        
        DO j = 1, ppm_dim
           
           IF ( j == 2 ) THEN
              
              xt(j,k) = xt(j,i) + len_phys(2)
              
           ELSE
              
              xt(j,k) = &
                   MODULO(xt(j,i) + shear_length(j,4), &
                   len_phys(j) )
              
           END IF
           
        END DO
        
     END IF
     
  END DO
  
  !----------------------------------------------------------
  ! north bc, if we are not using symmetry.
  ! With '='.
  !----------------------------------------------------------
  
  IF ( isymm == 0 ) THEN
     
     ymaxf = max_phys(2) 
     ymaxi = max_phys(2) - ghostsize
     
     DO i = 1, nghostplus
        
        IF ( xt(2,i) >  ymaxi .AND. &
             xt(2,i) <= ymaxf ) THEN
           
           k         = k + 1
           ighost(k) = ighost(i)
           
           DO j = 1, ppm_dim
              
              IF ( j == 2 ) THEN
                 
                 xt(j,k) = xt(j,i) - len_phys(2)
                 
              ELSE
                 
                 xt(j,k) = &
                      MODULO( xt(j,i) + shear_length(j,3),&
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
  
