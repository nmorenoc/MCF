  !----------------------------------------------------------
  ! (Re)allocate memory for the periodic ghosts
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
  ! periodic ghosts in the y-direction.
  !----------------------------------------------------------
  
  k  = nghostplus

  !----------------------------------------------------------
  ! first those at the south boundary.
  ! No '=' for inequality check at ghostize away.
  !----------------------------------------------------------
    
  yminf = min_phys(2) 
  ymini = min_phys(2) + ghostsize
  
  DO i = 1, nghostplus
     
     IF ( xt(2,i) >= yminf .AND. &
          xt(2,i) <  ymini ) THEN
        
        k         = k + 1
        ighost(k) = ighost(i)
        xt(1,k)   = xt(1,i) 
        xt(2,k)   = xt(2,i) + len_phys(2)
        
        IF ( ppm_dim == 3 ) THEN
           
           xt(3,k)   = xt(3,i)
           
        END IF
        
     END IF
     
  END DO
  
  IF ( isymm == 0 ) THEN
     
     !-------------------------------------------------------
     ! then north bc, but only if we are not using symmetry.
     ! No  '=' for inequality check at ghostzie away.
     !-------------------------------------------------------
     
     ymaxf = max_phys(2) 
     ymaxi = max_phys(2) - ghostsize
     
     DO i = 1, nghostplus
        
        IF ( xt(2,i) >  ymaxi .AND. &
             xt(2,i) <= ymaxf ) THEN
           
           k         = k + 1
           ighost(k) = ighost(i)
           xt(1,k)   = xt(1,i)
           xt(2,k)   = xt(2,i) - len_phys(2)
           
           IF ( ppm_dim == 3 ) THEN
              
              xt(3,k)   = xt(3,i)
              
           END IF
           
        END IF
        
     END DO
     
  END IF
  
  !----------------------------------------------------------
  ! update the ghost counter.
  !----------------------------------------------------------
  
  nghostplus = k
  
  
