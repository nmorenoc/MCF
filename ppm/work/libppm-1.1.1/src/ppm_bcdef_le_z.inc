  !----------------------------------------------------------
  ! (Re)allocate memory for the Lees-Edwards ghosts.
  !---------------------------------------------------------
  
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
  ! Lees-Edwards ghosts in the z-direction.
  !----------------------------------------------------------
  
  k = nghostplus 
  
  !----------------------------------------------------------
  ! those at the down boundary.
  ! with '='.  
  !----------------------------------------------------------
  
  zminf = min_phys(3) 
  zmini = min_phys(3) + ghostsize
  
  DO i = 1, nghostplus
     
     IF ( xt(3,i) >= zminf .AND. &
          xt(3,i) < zmini ) THEN
        
        k         = k + 1
        ighost(k) = ighost(i)
        
        DO j = 1, ppm_dim
           
           IF ( j == 3 ) THEN
              
              xt(j,k) = xt(j,i) + len_phys(3)
              
           ELSE
              
              xt(j,k) = &
                   MODULO(xt(j,i) + shear_length(j,6), &
                   len_phys(j) )
              
           END IF
           
        END DO
        
     END IF
     
  END DO
  
  !----------------------------------------------------------
  ! up bc, if we are not using symmetry.
  ! with '='.
  !----------------------------------------------------------

  IF ( isymm == 0 ) THEN
     
     zmaxf = max_phys(3) 
     zmaxi = max_phys(3) - ghostsize
     
     DO i = 1, nghostplus
        
        IF ( xt(3,i) >  zmaxi .AND. &
             xt(3,i) <= zmaxf ) THEN
           
           k         = k + 1
           ighost(k) = ighost(i)
           
           DO j = 1, ppm_dim
              
              IF ( j == 3 ) THEN
                 
                 xt(j,k) = xt(j,i) - len_phys(3)
                 
              ELSE
                 
                 xt(j,k) = &
                      MODULO( xt(j,i) + shear_length(j,5), &
                      len_phys(j) )
                 
              END IF
              
           END DO
           
        END IF
        
     END DO
     
  END IF
  
  !----------------------------------------------------------
  ! update the ghost counter.
  !----------------------------------------------------------
  
  nghostplus = k
  
  
  
