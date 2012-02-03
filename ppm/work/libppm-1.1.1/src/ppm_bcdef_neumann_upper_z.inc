      !-------------------------------------------------------------------
      !  (Re)allocate memory for the periodic ghosts
      !-------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow_preserve
      ldu(1) = ppm_dim
      ldu(2) = 2*nghostplus
      CALL ppm_alloc(xt,ldu,iopt,info) 
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &             'xt',__LINE__,info)
         GOTO 9999
      ENDIF

      ldu(1) = ldu(2)
      CALL ppm_alloc(ighost,ldu,iopt,info) 
      IF (info.NE.0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',     &
     &              'ighost',__LINE__,info)
         GOTO 9999
      ENDIF
      IF (isymm.EQ.0) THEN
         !----------------------------------------------------------------
         !  then the east bc, but only if we are not using symmetry
         !----------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION
         zmaxf = ppm_max_physs(3,topoid) 
         zmaxi = ppm_max_physs(3,topoid) - ghostsize
#else
         zmaxf = ppm_max_physd(3,topoid) 
         zmaxi = ppm_max_physd(3,topoid) - ghostsize
#endif
         DO i=1,nghostplus
            IF  (xt(3,i).GT.zmaxi.AND.xt(3,i).LT.zmaxf) THEN
               k         = k + 1
               ighost(k) = ighost(i)
               !xt(1,k)   = xt(1,i) - len_phys(1)
               xt(1,k)   = xt(1,i) ! free
               xt(2,k)   = xt(2,i)
               xt(3,k)   = xt(3,i) + ghostsize
            ENDIF
         ENDDO
      ENDIF

      !-------------------------------------------------------------------
      !  update the ghost counter
      !-------------------------------------------------------------------
      nghostplus = k
