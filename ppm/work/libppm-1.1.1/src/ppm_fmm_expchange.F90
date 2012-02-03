#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_fmm_expchange
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine exchanges the expansion coefficents,
      !                 the radius and the centerofboxes between all processores
      !                 by all to all communication    
      !                
      !
      !  Input        : 
      !                 prec         (F) : dummy to determine precision
      !                 lda          (I) leading dimension of vector case
      !                                    
      !  Output       : info         (I) : return status, 0 on success
      !
      !  Remarks      : This routine has 3 separate sendrecv, ie. the data is
      !                        not packed before  sending
      !                       (no pack data -send - unpack data ) 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fmm_expchange.f,v $
      !  Revision 1.14  2006/09/04 18:34:46  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.12  2006/06/29 10:28:35  pchatela
      !  Added vector strengths support
      !
      !  Revision 1.11  2006/06/16 07:52:21  hiebers
      !  Added a new list of topo IDs (topoidlist) to prevent overwriting user defined
      !  topologies
      !
      !  Revision 1.10  2005/09/19 13:03:28  polasekb
      !  code cosmetics
      !
      !  Revision 1.9  2005/09/11 18:05:30  polasekb
      !  (final?) corrected version
      !  (also works parallel :-)
      !
      !  Revision 1.8  2005/09/11 11:44:20  polasekb
      !  also communicating radius and centerofbox
      !
      !  Revision 1.7  2005/08/30 08:48:16  polasekb
      !  removed debug output
      !
      !  Revision 1.6  2005/08/25 14:16:02  polasekb
      !  corrected size of send/recv buffers
      !  exchanged sendrank/recvrank
      !
      !  Revision 1.4  2005/08/23 14:30:28  polasekb
      !  now making difference between single/double precision
      !
      !  Revision 1.3  2005/08/23 07:56:24  polasekb
      !  added #ifdef USE_MPI where needed
      !
      !  Revision 1.2  2005/08/23 07:49:26  polasekb
      !  corrected error output
      !
      !  Revision 1.1  2005/05/27 07:57:54  polasekb
      !  initial implementation
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   (__KIND == __SINGLE_PRECISION && __DIM == __SFIELD)
      SUBROUTINE ppm_fmm_expchange_s_sf(prec,info)
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __SFIELD)
      SUBROUTINE ppm_fmm_expchange_d_sf(prec,info)
#elif (__KIND == __SINGLE_PRECISION && __DIM == __VFIELD)
      SUBROUTINE ppm_fmm_expchange_s_vf(lda,prec,info)
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __VFIELD)
      SUBROUTINE ppm_fmm_expchange_d_vf(lda,prec,info)
#endif


      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_fmm
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
     
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
#ifdef HAVE_MPI
      INCLUDE 'mpif.h'
#else
#include "fakempi.h"
#endif

      !-------------------------------------------------------------------------
      !  Precision
      !------------------------------------------------------------------------- 
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if   __DIM == __VFIELD
      INTEGER                 , INTENT(IN   ) :: lda
#endif      
      REAL(MK)                , INTENT(IN   ) :: prec
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! auxiliary variables
      INTEGER                              :: i,j,k,isub,m,n,l
      INTEGER                              :: nsend,nrecv,curtopoid,box
      INTEGER                              :: sendrank,recvrank
      INTEGER                              :: nsendexp,nrecvexp
      INTEGER                              :: nsendrad,nrecvrad
      INTEGER                              :: nsendcen,nrecvcen
      INTEGER                              :: iopt,level,cnt, topoid
      INTEGER                              :: tag2,istat
      INTEGER, DIMENSION(1)                :: ldu1
      INTEGER, DIMENSION(2)                :: ldu2
      INTEGER, DIMENSION(3)                :: ldu3
      INTEGER, DIMENSION(4)                :: ldu4      
      REAL(MK)                             :: t0
      
      ! parallelisation
#ifdef USE_MPI
      INTEGER, DIMENSION(MPI_STATUS_SIZE)  :: status
#endif
      
      ! fmm      
      REAL(MK),DIMENSION(:  ),     POINTER :: radius,recvrad,sendrad
      REAL(MK),DIMENSION(:,:),     POINTER :: recvcen,sendcen,centerofbox 

#if   __DIM == __SFIELD       
      COMPLEX(MK),DIMENSION(:,:,:),POINTER :: expansion,recvexp,sendexp
#else
      COMPLEX(MK),DIMENSION(:,:,:,:),POINTER :: expansion,recvexp,sendexp
#endif

      
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_fmm_expchange',t0,info)
      
      !-------------------------------------------------------------------------
      !  pointing to correct variables (single/double)
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
#if   __DIM == __SFIELD
      expansion   => expansion_s_sf
#else
      expansion   => expansion_s_vf     
#endif
      centerofbox => centerofbox_s
      radius      => radius_s
#else
#if   __DIM == __SFIELD
      expansion   => expansion_d_sf
#else
      expansion   => expansion_d_vf     
#endif
      centerofbox => centerofbox_d
      radius      => radius_d
#endif


      !-------------------------------------------------------------------------
      !  Allocate memory for the sendlist
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ppm_nsendlist = ppm_nproc
      ppm_nrecvlist = ppm_nproc
      
      ldu1(1)        = ppm_nsendlist
      CALL ppm_alloc(ppm_isendlist,ldu1,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange',     &
          &    'send list PPM_ISENDLIST',__LINE__,info)
          GOTO 9999
      ENDIF
      
      CALL ppm_alloc(ppm_irecvlist,ldu1,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange',     &
          &    'receive list PPM_IRECVLIST',__LINE__,info)
          GOTO 9999
      ENDIF

      tag2   = 200

      
      !-------------------------------------------------------------------------
      !  compute top level topology
      !-------------------------------------------------------------------------

#ifdef  __VECTOR
      DO i=1,nlevel
        IF (nbpl(i) .GE. ppm_nproc) THEN
           level = i
        ENDIF
      ENDDO
#else
      DO i=1,nlevel
        IF (nbpl(i) .GE. ppm_nproc) THEN
           level = i
           EXIT
        ENDIF
      ENDDO
#endif

      
      !-------------------------------------------------------------------------
      !   Compute nsendexp,nsendcen,nsendrad (nr of exp.,cen.,rad. to be sent)
      !-------------------------------------------------------------------------
      nsendexp = 0
      DO i=level,nlevel
         topoid = topoidlist(i)
         !--------------------------------------------------------------------
         !  Get the ppm internal topoid
         !--------------------------------------------------------------------
         curtopoid = ppm_internal_topoid(topoid)
         nsendexp = nsendexp + ppm_nsublist(curtopoid)
         nsendcen = nsendexp
         nsendrad = nsendexp
      ENDDO
      !-------------------------------------------------------------------------
      !   Set up own lists for sending to other processors
      !-------------------------------------------------------------------------
      ! expansion
#if   __DIM == __SFIELD      
      ldu3(1) = nsendexp
      ldu3(2) = order+1
      ldu3(3) = 2*order+1
      CALL ppm_alloc(sendexp,ldu3,iopt,info)
#else
      ldu4(1) = lda
      ldu4(2) = nsendexp
      ldu4(3) = order+1
      ldu4(4) = 2*order+1
      CALL ppm_alloc(sendexp,ldu4,iopt,info)
#endif 
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
      &       'error allocating sendexp',__LINE__,info)
      GOTO 9999
      ENDIF 

      ! centerofbox
      ldu2(1) = 3
      ldu2(2) = nsendcen
      CALL ppm_alloc(sendcen,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
      &       'error allocating sendcen',__LINE__,info)
      GOTO 9999
      ENDIF 

      ! radius
      ldu1(1) = nsendrad
      CALL ppm_alloc(sendrad,ldu1,iopt,info)
      
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
      &       'error allocating sendrad',__LINE__,info)
      GOTO 9999
      ENDIF 
      
      
      !-------------------------------------------------------------------------
      !  loop over levels of tree
      !-------------------------------------------------------------------------
      cnt = 0
      DO i=level,nlevel
         topoid = topoidlist(i)
         !----------------------------------------------------------------------
         !  Get the ppm internal topoid
         !--------------------------------------------------------------------
         curtopoid = ppm_internal_topoid(topoid)
         !----------------------------------------------------------------------
         !   loop over local subs and store in sendexp,sendcen,sendrad
         !----------------------------------------------------------------------
         
	 DO j=1,ppm_nsublist(curtopoid)
            
	    box = ppm_boxid(ppm_isublist(j,curtopoid),i)
            cnt = cnt + 1

#if         __DIM == __SFIELD
	    DO n=1,ldu3(2)
               DO m=1,ldu3(3)
                  sendexp(cnt,n,m) = expansion(box,n,m)       
               ENDDO
            ENDDO
#else
            DO n=1,ldu4(2)
               DO m=1,ldu4(3)
                  DO l=1,ldu4(1)
                     sendexp(l,cnt,n,m) = expansion(l,box,n,m)       
                  ENDDO
               ENDDO
            ENDDO
#endif
            DO n=1,3
               sendcen(n,cnt)   = centerofbox(n,box)
            ENDDO
            sendrad(cnt)     = radius(box)
         ENDDO
      ENDDO
      
      
      !-------------------------------------------------------------------------
      !   Initialize sendrank, recvrank, ppm_nsendlist, ppm_nrecvlist
      !-------------------------------------------------------------------------
      sendrank           = ppm_rank - 1
      recvrank           = ppm_rank + 1
      ppm_nsendlist      = 0
      ppm_nrecvlist      = 0

      !-------------------------------------------------------------------------
      !  Since we skip the local processor entirely, increment the pointers once
      !-------------------------------------------------------------------------
      sendrank                     = sendrank + 1
      recvrank                     = recvrank - 1
      ppm_nsendlist                = ppm_nsendlist + 1
      ppm_isendlist(ppm_nsendlist) = sendrank
      ppm_nrecvlist                = ppm_nrecvlist + 1
      ppm_irecvlist(ppm_nrecvlist) = recvrank

      !-------------------------------------------------------------------------
      !  Loop over all processors but skip the processor itself
      !-------------------------------------------------------------------------
      DO i=2,ppm_nproc
	 !----------------------------------------------------------------------
         !  compute the next processor
         !----------------------------------------------------------------------
	 sendrank = sendrank + 1
         IF (sendrank.GT.ppm_nproc-1) sendrank = sendrank - ppm_nproc 
         recvrank = recvrank - 1
         IF (recvrank.LT.          0) recvrank = recvrank + ppm_nproc 

	 !----------------------------------------------------------------------
         !  Store the processor to which we will send to
         !----------------------------------------------------------------------
	 ppm_nsendlist                = ppm_nsendlist + 1 
         ppm_isendlist(ppm_nsendlist) = sendrank

         !----------------------------------------------------------------------
         !  Store the processor to which we will recv from
         !----------------------------------------------------------------------
	 ppm_nrecvlist                = ppm_nrecvlist + 1 
         ppm_irecvlist(ppm_nrecvlist) = recvrank
  
         
	 !----------------------------------------------------------------------
         !  reset counter for nr of exp.coeff.,rad,centers to be received
         !----------------------------------------------------------------------
	 nrecvexp                     = 0
         nrecvcen                     = 0
         nrecvrad                     = 0
         
	 !----------------------------------------------------------------------
         !  loop over all topologies and check all subs 
         !----------------------------------------------------------------------
         DO j=level,nlevel
            topoid = topoidlist(j)
            !-------------------------------------------------------------------
            !  Get the ppm internal topoid
            !------------------------------------------------------------------
            curtopoid = ppm_internal_topoid(topoid)
            DO isub=1,ppm_nsubs(curtopoid)
              !-----------------------------------------------------------------
              !  Check if they belong to the processor from where we will 
              !  receive data
              !-----------------------------------------------------------------
              IF (ppm_subs2proc(isub,curtopoid) .EQ. recvrank) THEN
                 !--------------------------------------------------------------
                 !  If yes, increase counter for correct allocation
                 !--------------------------------------------------------------
                 nrecvexp = nrecvexp + 1
                 nrecvcen = nrecvcen + 1
                 nrecvrad = nrecvrad + 1
              ELSE
                 !--------------------------------------------------------------
                 !  will be exchanged in another round
                 !--------------------------------------------------------------
              ENDIF 
           ENDDO
         ENDDO
         
	 
	 !----------------------------------------------------------------------
         !  Allocate our recv-arrays, exp.,rad., centers to be received
         !----------------------------------------------------------------------
         iopt   = ppm_param_alloc_fit

	 ! expansion
#if      __DIM == __SFIELD
         ldu3(1) = nrecvexp
         ldu3(2) = order+1
         ldu3(3) = 2*order+1
         CALL ppm_alloc(recvexp,ldu3,iopt,info)
#else
	 ldu4(1) = lda
         ldu4(2) = nrecvexp
         ldu4(3) = order+1
         ldu4(4) = 2*order+1
         CALL ppm_alloc(recvexp,ldu4,iopt,info)
#endif

         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
         &       'error allocating recvexp',__LINE__,info)
         GOTO 9999
         ENDIF

	 ! centerofbox
         ldu2(1) = 3
         ldu2(2) = nrecvcen
         CALL ppm_alloc(recvcen,ldu2,iopt,info)

         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
         &       'error allocating recvcen',__LINE__,info)
         GOTO 9999
         ENDIF

	 ! radius
         ldu1(1) = nrecvrad
         CALL ppm_alloc(recvrad,ldu1,iopt,info)
 
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fmm_expchange', &
         &       'error allocating recvrad',__LINE__,info)
         GOTO 9999
         ENDIF

	 
         !----------------------------------------------------------------------
         !  receive the expansions
         !----------------------------------------------------------------------
	 nsend = nsendexp*(order+1)*(2*order+1)
         nrecv = nrecvexp*(order+1)*(2*order+1)
#ifdef   USE_MPI
#if      __KIND == __SINGLE_PRECISION
         CALL MPI_SendRecv(sendexp,nsend,MPI_COMPLEX,sendrank,tag2, &
         &                 recvexp,nrecv,MPI_COMPLEX,recvrank,tag2, &
         &                 ppm_comm,status,info)
#else
         CALL MPI_SendRecv(sendexp,nsend,MPI_DOUBLE_COMPLEX,sendrank,tag2,&
         &                 recvexp,nrecv,MPI_DOUBLE_COMPLEX,recvrank,tag2,&
         &                 ppm_comm,status,info)
#endif
#endif
         !----------------------------------------------------------------------
         !  receive the centers
         !----------------------------------------------------------------------
	 nsend = nsendcen*3
         nrecv = nrecvcen*3
#ifdef   USE_MPI
         CALL MPI_SendRecv(sendcen,nsend,ppm_mpi_kind,sendrank,tag2, &
         &                 recvcen,nrecv,ppm_mpi_kind,recvrank,tag2, &
         &                 ppm_comm,status,info)
#endif
         !----------------------------------------------------------------------
         !  receive the radius
         !---------------------------------------------------------------------- 
	 nsend = nsendrad
         nrecv = nrecvrad
#ifdef USE_MPI
         CALL MPI_SendRecv(sendrad,nsend,ppm_mpi_kind,sendrank,tag2,&
         &                 recvrad,nrecv,ppm_mpi_kind,recvrank,tag2,&
         &                 ppm_comm,status,info)
#endif
         !----------------------------------------------------------------------
         !  store the received data
         !  loop over all topologies and check all subs
         !----------------------------------------------------------------------
         cnt = 0
         DO j=level,nlevel
            topoid = topoidlist(j)
            !-------------------------------------------------------------------
            !  Get the ppm internal topoid
            !-----------------------------------------------------------------
            curtopoid = ppm_internal_topoid(topoid)
            DO isub=1,ppm_nsubs(curtopoid)
              !-----------------------------------------------------------------
              !  Check if the sub belongs to the processor from where we 
              !  received data
              !-----------------------------------------------------------------
	      IF (ppm_subs2proc(isub,curtopoid) .EQ. recvrank) THEN
		 box = ppm_boxid(isub,j)
		 !--------------------------------------------------------------
                 !  If yes, store
                 !--------------------------------------------------------------
                 cnt = cnt + 1
#if              __DIM == __SFIELD
                 DO n=1,ldu3(2)
                    DO m=1,ldu3(3)
                       expansion(box,n,m) = recvexp(cnt,n,m)
                    ENDDO
                 ENDDO
#else
                 DO n=1,ldu4(2)
                    DO m=1,ldu4(3)
                       DO l=1,ldu4(1)
                          expansion(l,box,n,m) = recvexp(l,cnt,n,m)
                       ENDDO
                    ENDDO
                 ENDDO

#endif
                 DO n=1,3
                    centerofbox(n,box) = recvcen(n,cnt)
                 ENDDO
                 radius(box)        = recvrad(cnt)
              ENDIF
           ENDDO
         ENDDO
      ENDDO ! end loop over nproc

      
      !-------------------------------------------------------------------------
      !  Deallocate the memory for the lists
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      ldu3(1) = 0
      ldu3(2) = 0
      ldu3(3) = 0

      istat   = 0
      CALL ppm_alloc(recvexp,ldu3,iopt,info)
      istat = istat + info
      CALL ppm_alloc(recvcen,ldu3,iopt,info)
      istat = istat + info
      CALL ppm_alloc(recvrad,ldu3,iopt,info)
      istat = istat + info
      IF (istat.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_fmm_expchange',     &
         &    'recvexp',__LINE__,info)
      ENDIF

      istat   = 0
      CALL ppm_alloc(sendexp,ldu3,iopt,info)
      istat = istat + info
      CALL ppm_alloc(sendcen,ldu3,iopt,info)
      istat = istat + info
      CALL ppm_alloc(sendrad,ldu3,iopt,info)
      istat = istat + info
      IF (istat.NE.0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_fmm_expchange',     &
         &    'sendexp',__LINE__,info)
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fmm_expchange',t0,info)
      RETURN
#if   (__KIND == __SINGLE_PRECISION && __DIM == __SFIELD)
      END SUBROUTINE ppm_fmm_expchange_s_sf
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __SFIELD)
      END SUBROUTINE ppm_fmm_expchange_d_sf
#elif (__KIND == __SINGLE_PRECISION && __DIM == __VFIELD)
      END SUBROUTINE ppm_fmm_expchange_s_vf
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __VFIELD)
      END SUBROUTINE ppm_fmm_expchange_d_vf
#endif
