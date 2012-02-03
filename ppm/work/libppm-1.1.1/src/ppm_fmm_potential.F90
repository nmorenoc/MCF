#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_fmm_potential
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Compute the potential at the target points using
      !                 the expansion coefficients
      !                 In the parallel case:
      !                 Calls ppm_fmm_pretraverse and ppm_fmm_expchange
      !                 Maps the target points onto the leaf topolgy
      !                 
      !  Input        : xpunord(:,:) (F) the position of the field points
      !                 wpunord(:)   (F) the strength of the field points
      !                 tp(:,:)      (F) the target points
      !                 theta        (F) acceptance factor
      !
      !  Input/output :     
      !                 Np           (I) the number of field points.
      !                 Ntp          (I) the number of target points
      !
      !  Output       : potential(:) (F) the multipole expansion potential
      !                                  for each point
      !                                  size 1:Np
      !                 info         (I) return status. 0 upon success.
      !
      !  Remarks      :  
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fmm_potential.f,v $
      !  Revision 1.30  2006/09/04 18:34:46  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.29  2006/06/29 10:28:35  pchatela
      !  Added vector strengths support
      !
      !  Revision 1.28  2006/06/20 15:15:51  hiebers
      !  adjusted arguments in call of ppm_fmm_pretraverse
      !
      !  Revision 1.27  2006/06/16 07:52:21  hiebers
      !  Added a new list of topo IDs (topoidlist) to prevent overwriting user defined
      !  topologies
      !
      !  Revision 1.26  2006/06/09 07:56:46  pchatela
      !  Bugfix: weights wp were integers!
      !
      !  Revision 1.25  2005/09/19 13:03:29  polasekb
      !  code cosmetics
      !
      !  Revision 1.24  2005/09/12 09:14:10  hiebers
      !  added mapping of target points
      !
      !  Revision 1.23  2005/09/11 18:05:31  polasekb
      !  (final?) corrected version
      !  (also works parallel :-)
      !
      !  Revision 1.22  2005/09/10 07:50:04  polasekb
      !  changed init of stack for parallel version
      !
      !  Revision 1.21  2005/09/05 06:23:57  polasekb
      !  corrected variable initialisation
      !
      !  Revision 1.20  2005/08/29 15:18:00  polasekb
      !  bugfix when computing direct way
      !
      !  Revision 1.19  2005/08/25 13:52:10  polasekb
      !  mapping corrected from wp to wpunord
      !
      !  Revision 1.18  2005/08/23 14:35:20  polasekb
      !  changed call to ppm_fmm_pretraverse
      !
      !  Revision 1.17  2005/08/23 14:30:05  polasekb
      !  changed call to ppm_fmm_expchange
      !
      !  Revision 1.16  2005/08/23 14:16:50  polasekb
      !  corrected wpunord and wp
      !
      !  Revision 1.15  2005/08/11 15:12:04  polasekb
      !  fixed indices of Outer
      !
      !  Revision 1.14  2005/08/11 13:32:34  polasekb
      !  added variable theta
      !
      !  Revision 1.13  2005/08/08 13:35:52  polasekb
      !  deallocate some local variables
      !
      !  Revision 1.12  2005/08/04 16:04:51  polasekb
      !  removed some data allocation
      !
      !  Revision 1.11  2005/07/29 14:06:52  polasekb
      !  changed check of eqalness to ppm_myeps
      !
      !  Revision 1.10  2005/07/29 12:36:13  polasekb
      !  changed diagonal to radius
      !
      !  Revision 1.9  2005/07/27 14:59:24  polasekb
      !  now using constants from data file
      !
      !  Revision 1.8  2005/07/25 15:01:20  polasekb
      !  bugfix with indices
      !
      !  Revision 1.7  2005/07/25 14:40:31  polasekb
      !  adapted computation of the potential
      !
      !  Revision 1.6  2005/07/21 12:42:07  polasekb
      !  bugfix in allocating an array
      !
      !  Revision 1.5  2005/07/21 08:26:09  polasekb
      !  changed function call, now different target 
      !  points and field points can be
      !  specified by the user
      !
      !  Revision 1.4  2005/06/02 19:18:51  polasekb
      !  corrected syntax error
      !
      !  Revision 1.3  2005/05/30 09:37:01  polasekb
      !  bugfix: corrected call to ppm_util_cart2sph
      !
      !  Revision 1.2  2005/05/27 12:42:48  polasekb
      !  initialized further arrays
      !
      !  Revision 1.1  2005/05/27 08:01:23  polasekb
      !  initial implementation
      !  TODO: remove debug output
      !
      !  Revision 0  2005/01/16 15:59:14 polasekb
      !  start
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   (__KIND == __SINGLE_PRECISION && __DIM == __SFIELD)
      SUBROUTINE ppm_fmm_potential_s_sf(xpunord,wpunord,Np,tp,Ntp,theta, &
                 &                      potential,info)
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __SFIELD)
      SUBROUTINE ppm_fmm_potential_d_sf(xpunord,wpunord,Np,tp,Ntp,theta, &
                 &                      potential,info)
#elif (__KIND == __SINGLE_PRECISION && __DIM == __VFIELD)
      SUBROUTINE ppm_fmm_potential_s_vf(xpunord,wpunord,lda,Np,tp,Ntp,theta,&
                 &                      potential,info)
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __VFIELD)
      SUBROUTINE ppm_fmm_potential_d_vf(xpunord,wpunord,lda,Np,tp,Ntp,theta,&
                 &                      potential,info)
#endif


      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      
      USE ppm_module_data
      USE ppm_module_data_fmm
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_fmm_expchange
      USE ppm_module_fmm_pretraverse
      USE ppm_module_map
      USE ppm_module_map_part_get_sub
      USE ppm_module_mktopo
      USE ppm_module_substart
      USE ppm_module_substop 
      USE ppm_module_topo
      USE ppm_module_topo_box2subs
      USE ppm_module_util_cart2sph
      USE ppm_module_write
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
      REAL(MK), DIMENSION(:,:), POINTER       :: xpunord
      INTEGER                 , INTENT(INOUT) :: Np
      REAL(MK), DIMENSION(:,:), POINTER       :: tp
      INTEGER                 , INTENT(INOUT) :: Ntp
      REAL(MK)                , INTENT(IN   ) :: theta
      INTEGER                 , INTENT(  OUT) :: info

#if   __DIM == __SFIELD
      REAL(MK), DIMENSION(:), POINTER         :: wpunord
      REAL(MK), DIMENSION(:), POINTER         :: potential
#else
      REAL(MK), DIMENSION(:,:), POINTER         :: wpunord
      INTEGER                                   :: lda
      REAL(MK), DIMENSION(:,:), POINTER         :: potential
#endif
     
      
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      LOGICAL                              :: check,drct,OK
      INTEGER                              :: i,j,k,l,cnt,iopt,m,n
      INTEGER                              :: pcount,ccount
      INTEGER                              :: mapt,Mpart,root,istat
      INTEGER                              :: first,last,level
      INTEGER                              :: stackpointer,curbox
      INTEGER                              :: pexp,isymm
      INTEGER ,DIMENSION(1)                :: ldu1
      INTEGER ,DIMENSION(2)                :: ldu2 
      INTEGER ,DIMENSION(:  ), POINTER     :: newlpdx,stack
      REAL(MK)                             :: thetap,eps,angle,reci 
      REAL(MK)                             :: sine,cosine,val,prod 
      REAL(MK),DIMENSION(1)                :: curboxrho,curboxphi,curboxtheta
      REAL(MK),DIMENSION(:,:),     POINTER :: min_box,max_box
      REAL(MK),DIMENSION(:,:),     POINTER :: min_sub,max_sub
      COMPLEX(MK),PARAMETER                :: CI=(0.0_MK,1.0_MK)
      CHARACTER(LEN=ppm_char)              :: cbuf
      REAL(MK)                             :: dx,dy,dz,dist,rad
      ! parallelisation
      REAL(MK)                             :: t0,ghostsize,cutoff
      INTEGER                              :: topoid
      INTEGER ,DIMENSION(:  ), POINTER     :: part_subtop
            
      ! fmm 
      REAL(MK),DIMENSION(:  ),     POINTER :: fracfac,boxcost
      REAL(MK),DIMENSION(:,:),     POINTER :: sqrtfac,xp,Anm
      REAL(MK),DIMENSION(:  ), POINTER     :: radius
      REAL(MK),DIMENSION(:,:), POINTER     :: Pnm,centerofbox 
      COMPLEX(MK),DIMENSION(:,:),  POINTER :: Ynm 
      COMPLEX(MK),DIMENSION(:,:),  POINTER :: Outer             

#if   __DIM == __SFIELD
      REAL ,DIMENSION(:  ), POINTER        :: wp      
      COMPLEX(MK),DIMENSION(:,:,:),POINTER :: expansion
#else 
      REAL ,DIMENSION(:,:), POINTER          :: wp
      COMPLEX(MK),DIMENSION(:,:,:,:),POINTER :: expansion           
#endif     
      
      
      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------
      CALL substart('ppm_fmm_potential',t0,info)
      
      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN  
            IF (Np .LT. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_fmm_potential',   &
      &               'number of particles must be > 0 !',__LINE__,info)
               GOTO 9999
            ENDIF
            IF (.NOT. ppm_fmm_initialized) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_fmm_potential',   &
      &               'Please call ppm_fmm_init first',__LINE__,info)
               GOTO 9999
            ENDIF

      ENDIF

      
      !-------------------------------------------------------------------------
      ! Check precision and pointing tree data to correct variables
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      min_box      => min_box_s
      max_box      => max_box_s
      boxcost      => boxcost_s
      centerofbox  => centerofbox_s
      radius       => radius_s
#if   __DIM == __SFIELD
      expansion    => expansion_s_sf
#else
      expansion    => expansion_s_vf
#endif      
      sqrtfac      => sqrtfac_s
      fracfac      => fracfac_s
      Anm          => Anm_s
      eps          = ppm_myepss
      Ynm          => Ynm_s
      Pnm          => Pnm_s
      Outer        => Outer_s
#else
      min_box      => min_box_d
      max_box      => max_box_d
      boxcost      => boxcost_d
      centerofbox  => centerofbox_d
      radius       => radius_d
#if   __DIM == __SFIELD
      expansion    => expansion_d_sf
#else
      expansion    => expansion_d_vf
#endif
      sqrtfac      => sqrtfac_d
      fracfac      => fracfac_d
      Anm          => Anm_d
      eps          = ppm_myepsd
      Ynm          => Ynm_d
      Pnm          => Pnm_d
      Outer        => Outer_d
#endif

      IF (ppm_nproc .GT. 1) THEN
      !choose lowest level of the tree as topology
      topoid = topoidlist(nlevel)
      

      !-------------------------------------------------------------------------
      ! Call fmm_expchange to communicate all expansions
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
      CALL ppm_fmm_expchange(t0,info)
#else      
      CALL ppm_fmm_expchange(lda,t0,info)
#endif
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to call expchange.',info)
          GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      ! Map target points
      !-------------------------------------------------------------------------
      mapt = ppm_param_map_global
      CALL ppm_map_part(tp,ppm_dim,Ntp,Mpart,topoid,mapt,info)   ! positions
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to start global mapping.',info)
          GOTO 9999
      ENDIF
      
      mapt = ppm_param_map_send
      CALL ppm_map_part(tp,ppm_dim,Ntp,Mpart,topoid,mapt,info)   ! positions
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to start global mapping.',info)
          GOTO 9999
      ENDIF
      
      mapt = ppm_param_map_pop
      CALL ppm_map_part(tp,ppm_dim,Ntp,Mpart,topoid,mapt,info)   ! positions
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to start global mapping.',info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Store new number of particles
      !-------------------------------------------------------------------------
      Ntp = Mpart
      
      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &          'Done mapping target points.',info)
          WRITE(cbuf,'(A,I8)') 'Local number of target points now:',Ntp
          CALL ppm_write(ppm_rank,'ppm_fmm_potential',cbuf,info)
      ENDIF

      
      !-------------------------------------------------------------------------
      !  Check that particles have been mapped correctly
      !-------------------------------------------------------------------------

      IF(ppm_debug.GT.0)THEN
         CALL ppm_topo_check(tp,Ntp,OK,info)
         IF (info .NE. 0) THEN
            CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
         &    'Failed to check topology.',info)
         ENDIF
      
         IF (.NOT.OK) THEN
            CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'target points not mapped correctly!',info)
            GOTO 9999
         ENDIF
      ENDIF
      
      !-------------------------------------------------------------------------
      ! Call pretraversal routine to build communication lists
      !-------------------------------------------------------------------------
      CALL ppm_fmm_pretraverse(tp,Ntp,nlevel,theta,ccount,&
      &                       part_subtop,pcount,info)
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to call pretraverse.',info)
          GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      ! Get particles on local processor stored in part_subtop
      !-------------------------------------------------------------------------
      CALL ppm_map_part_get_sub(part_subtop,pcount,topoid,xpunord,Np,info)
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to call map_part_get_sub.',info)
          GOTO 9999
      ENDIF

      
      !-------------------------------------------------------------------------
      ! Push the weights and boxpart 
      !-------------------------------------------------------------------------
      isymm  = 0
      cutoff = 1.0_MK ! can be any number > 0

      mapt = ppm_param_map_push
      
#if   __DIM == __SFIELD
      CALL ppm_map_part_ghost(wpunord,Np,Mpart,isymm,cutoff,mapt,info) !strengths
#else      
      CALL ppm_map_part_ghost(wpunord,lda,Np,Mpart,isymm,cutoff,mapt,info) !strengths
#endif      
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF
      
      CALL ppm_map_part_ghost(boxpart,Np,Mpart,isymm,cutoff,mapt,info) !boxpart
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF
            
      mapt = ppm_param_map_send
      
      CALL ppm_map_part_ghost(boxpart,Np,Mpart,isymm,cutoff,mapt,info)   ! send
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to send particles.',info)
          GOTO 9999
      ENDIF
      
      mapt = ppm_param_map_pop
      
      CALL ppm_map_part_ghost(boxpart,Np,Mpart,isymm,cutoff,mapt,info)  !boxpart
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF      
      
#if   __DIM == __SFIELD      
      CALL ppm_map_part_ghost(wpunord,Np,Mpart,isymm,cutoff,mapt,info) !strengths
#else
      CALL ppm_map_part_ghost(wpunord,lda,Np,Mpart,isymm,cutoff,mapt,info) !strengths
#endif      
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to pop strengths.',info)
          GOTO 9999
      ENDIF !positions
      
      CALL ppm_map_part_ghost(xpunord,ppm_dim,Np,Mpart,isymm,cutoff,mapt,info)
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &    'Failed to pop positions.',info)
          GOTO 9999
      ENDIF

      IF (ppm_debug .GT. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
          &          'Done mapping ghost particles.',info)
          WRITE(cbuf,'(A,I8)') 'Received ghost particles:',Mpart-Np
          CALL ppm_write(ppm_rank,'ppm_fmm_potential',cbuf,info)
      ENDIF

      
      !-------------------------------------------------------------------------
      ! Sort new (ghost) particles 
      !-------------------------------------------------------------------------
      ldu1(1) = Mpart
      CALL ppm_alloc(newlpdx,ldu1,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_potential', &
      &       'error allocating newlpdx',__LINE__,info)
      GOTO 9999
      ENDIF

      newlpdx(1:Np) = lpdx(1:Np)

      cnt = Np
      DO i=1,nbox
         check = .TRUE.
         DO j=Np+1,Mpart
            IF (boxpart(j) .EQ. i) THEN
               cnt = cnt + 1
               IF (((lhbx(1,i) .EQ. 1) .OR. (lhbx(1,i) .EQ. Np+1)) .AND. &
                    (check)) THEN
                  check = .FALSE.
                  lhbx(1,i) = cnt
                  lhbx(2,i) = cnt
                  newlpdx(cnt) = j
               ELSE
                  lhbx(2,i) = lhbx(2,i) + 1
                  newlpdx(cnt) = j
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
      !-------------------------------------------------------------------------
      ! ATTENTION: now tree lists lhbx and lpdx are not valid anymore
      !-------------------------------------------------------------------------
      
      ENDIF !ppm_nproc .GT. 1

      
      !-------------------------------------------------------------------------
      ! Allocate array for potentials of target points
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit

#if   __DIM == __SFIELD
      ldu1(1) = Ntp
      CALL ppm_alloc(potential,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_potential', &
      &       'error allocating potential',__LINE__,info)
      GOTO 9999
      ENDIF
#else
      ldu2(1) = lda
      ldu2(2) = Ntp
      CALL ppm_alloc(potential,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_potential', &
      &       'error allocating potential',__LINE__,info)
      GOTO 9999
      ENDIF
#endif


      !-------------------------------------------------------------------------
      ! Initialize arrays
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
      DO i=1,Ntp
         potential(i) = 0.0_MK
      ENDDO
#else 
      DO i=1,Ntp
         DO j=1,lda
            potential(j,i) = 0.0_MK
         ENDDO
      ENDDO
#endif


      !-------------------------------------------------------------------------
      ! allocate and initialize further arrays
      !-------------------------------------------------------------------------
      istat = 0
      
      ldu1(1) = nbox
      CALL ppm_alloc(stack,ldu1,iopt,info)
      istat = istat + info
 
      ldu2(1) = 3
      IF (ppm_nproc .GT. 1) THEN
         ldu2(2) = Mpart
      ELSE
         ldu2(2) = Np
      ENDIF
      CALL ppm_alloc(xp,ldu2,iopt,info)
      istat = istat + info

#if   __DIM == __SFIELD
      IF (ppm_nproc .GT. 1) THEN
         ldu1(1) = Mpart
      ELSE
         ldu1(1) = Np
      ENDIF
      CALL ppm_alloc(wp,ldu1,iopt,info)
#else 
      ldu2(1) = lda     
      IF (ppm_nproc .GT. 1) THEN
         ldu2(2) = Mpart
      ELSE
         ldu2(2) = Np
      ENDIF
      CALL ppm_alloc(wp,ldu2,iopt,info)
#endif

      istat = istat + info
      
      IF (istat .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_potential', &
      &       'error allocating variables',__LINE__,info)
      GOTO 9999
      ENDIF

      stack = (0)
      xp    = (0.0_MK,0.0_MK)

#if   __DIM == __SFIELD
      DO i =1,size(wp,1)
         wp(i) = 0.0_MK
      ENDDO
#else
      DO i =1,size(wp,2)
         DO j =1,lda
            wp(i,j) = 0.0_MK
         ENDDO
      ENDDO
#endif
      

      !-------------------------------------------------------------------------
      ! Find the root of the tree (serial) and find the
      ! top level of tree (parallel)
      !-------------------------------------------------------------------------
#ifdef   __VECTOR
      IF (ppm_nproc .GT. 1) THEN        
         ! finding top level
         DO i=1,nlevel
           IF (nbpl(i) .GE. ppm_nproc) THEN
              level = i
           ENDIF
         ENDDO
      ELSE
        IF (parent(1) .EQ. ppm_param_undefined) THEN
           root = 1
        ELSE
           DO i=1,nbox
              IF (parent(i) .EQ. ppm_param_undefined) THEN
                 root = i
              ENDIF
           ENDDO
        ENDIF
      ENDIF

#else
      IF (ppm_nproc .GT. 1) THEN        
         ! finding top level
         DO i=1,nlevel
           IF (nbpl(i) .GE. ppm_nproc) THEN
              level = i
              EXIT
           ENDIF
         ENDDO
      ELSE
        IF (parent(1) .EQ. ppm_param_undefined) THEN
           root = 1
        ELSE
           DO i=1,nbox
              IF (parent(i) .EQ. ppm_param_undefined) THEN
                 root = i
                 EXIT
              ENDIF
           ENDDO
        ENDIF
      ENDIF
#endif

      !-------------------------------------------------------------------------
      ! order the particles according to the tree
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
      IF(ppm_dim.EQ.2)THEN
         IF (ppm_nproc .GT. 1) THEN
            DO i=1,nbox
               first = lhbx(1,i)
               last  = lhbx(2,i)
               DO j=first,last
                  xp(1,j) = xpunord(1,newlpdx(j))
                  xp(2,j) = xpunord(2,newlpdx(j))
                  wp(j)   = wpunord(newlpdx(j))
               ENDDO
            ENDDO
         ELSE
            DO i=1,nbox
               first = lhbx(1,i)
               last  = lhbx(2,i)
               DO j=first,last
                  xp(1,j) = xpunord(1,lpdx(j))
                  xp(2,j) = xpunord(2,lpdx(j))
                  wp(j)   = wpunord(lpdx(j))
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      IF(ppm_dim.EQ.3)THEN
         IF (ppm_nproc .GT. 1) THEN
            DO i=1,nbox
               first = lhbx(1,i)
               last  = lhbx(2,i)
               DO j=first,last
                  xp(1,j) = xpunord(1,newlpdx(j))
                  xp(2,j) = xpunord(2,newlpdx(j))
                  xp(3,j) = xpunord(3,newlpdx(j))
                  wp(j)   = wpunord(newlpdx(j))
               ENDDO
            ENDDO
         ELSE
            DO i=1,nbox
               first = lhbx(1,i)
               last  = lhbx(2,i)
               DO j=first,last
                  xp(1,j) = xpunord(1,lpdx(j))
                  xp(2,j) = xpunord(2,lpdx(j))
                  xp(3,j) = xpunord(3,lpdx(j))
                  wp(j)   = wpunord(lpdx(j))
               ENDDO
            ENDDO
         ENDIF
      ENDIF
#else 
      IF(ppm_dim.EQ.2)THEN
         IF (ppm_nproc .GT. 1) THEN
            DO i=1,nbox
               first = lhbx(1,i)
               last  = lhbx(2,i)
               DO j=first,last
                  xp(1,j) = xpunord(1,newlpdx(j))
                  xp(2,j) = xpunord(2,newlpdx(j))
                  DO l=1,lda
                     wp(l,j) = wpunord(l,newlpdx(j))
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO i=1,nbox
               first = lhbx(1,i)
               last  = lhbx(2,i)
               DO j=first,last
                  xp(1,j) = xpunord(1,lpdx(j))
                  xp(2,j) = xpunord(2,lpdx(j))
                  DO l=1,lda
                     wp(l,j) = wpunord(l,lpdx(j))
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      IF(ppm_dim.EQ.3)THEN
         IF (ppm_nproc .GT. 1) THEN
            DO i=1,nbox
               first = lhbx(1,i)
               last  = lhbx(2,i)
               DO j=first,last
                  xp(1,j) = xpunord(1,newlpdx(j))
                  xp(2,j) = xpunord(2,newlpdx(j))
                  xp(3,j) = xpunord(3,newlpdx(j))
                  DO l=1,lda
                     wp(l,j) = wpunord(l,newlpdx(j))
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            DO i=1,nbox
               first = lhbx(1,i)
               last  = lhbx(2,i)
               DO j=first,last
                  xp(1,j) = xpunord(1,lpdx(j))
                  xp(2,j) = xpunord(2,lpdx(j))
                  xp(3,j) = xpunord(3,lpdx(j))
                  DO l=1,lda
                     wp(l,j) = wpunord(l,lpdx(j))
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF

#endif


      !-------------------------------------------------------------------------
      ! Compute the potential for the target points
      ! Ntp: number of target points
      !-------------------------------------------------------------------------
      DO i=1,Ntp
         IF (ppm_nproc .GT. 1) THEN        
           
	   ! init stack parallel
           stackpointer = 1
           cnt = 0

           !--------------------------------------------------------------------
           ! Collect all boxes at the highest level (doesnt vectorize)
           !--------------------------------------------------------------------
	   DO j=1,nbox
             IF (blevel(j) .EQ. level) THEN
	       stack(stackpointer) = j 
               stackpointer = stackpointer + 1
               cnt = cnt +1
               IF (cnt .EQ. nbpl(level)) THEN
                  EXIT
               ENDIF
             ENDIF
           ENDDO
         ELSE
          
	  ! init stack serial
            stackpointer = 1
            stack(stackpointer) = root
            stackpointer = stackpointer + 1
         ENDIF
         
	 DO WHILE (stackpointer .GT. 1)
	   !pop top box
            stackpointer = stackpointer - 1
            curbox = stack(stackpointer)
           
           
            dx = tp(1,i) - centerofbox(1,curbox)
            dy = tp(2,i) - centerofbox(2,curbox)
            dz = tp(3,i) - centerofbox(3,curbox)
            dist = SQRT(dx*dx + dy*dy + dz*dz)
           
	   
	   !--------------------------------------------------------------------
           ! Checking Barnes-Hut Criterium
           !--------------------------------------------------------------------
	   drct = .FALSE.
           
           IF (radius(curbox) .LE. 0.0_MK) THEN
              !only one particle in box, do direct computation
              drct = .TRUE.
           ENDIF
           
	   IF ((dist/(2*radius(curbox)) .GT. theta) .AND. (.NOT. drct)) THEN
	     !-----------------------------------------------------------------
             !  far enough, compute part-box interaction
             !-----------------------------------------------------------------
	     IF (ppm_debug .GT. 0) THEN
                  CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
                  &    'far enough',info)
             ENDIF
             
	     
	     !------------------------------------------------------------------
             ! TODO:
             ! Computing the expansion order pexp according to Wang
             !------------------------------------------------------------------
	     !thetap = 0.75_MK+0.2_MK*(order-pexp)+0.05_MK*(order-pexp)**2
             pexp = order
             
	     !DO WHILE((thetap .LE. dist/radius(curbox)).AND.(pexp .GT. 3))
             !   pexp = pexp - 1
             !   thetap = 0.75_MK+0.2_MK*(order-pexp)+0.05_MK*(order-pexp)**2
             !ENDDO

             CALL ppm_util_cart2sph(tp(1,i:i),tp(2,i:i),tp(3,i:i),1, &
                 & centerofbox(1,curbox),centerofbox(2,curbox), &
                 & centerofbox(3,curbox), &
                 & curboxrho,curboxtheta,curboxphi,info)
             IF (info .NE. 0) THEN
                CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_potential', &
                     & 'Failed calling util_cart2sph',__LINE__,info)
             ENDIF
             
	     !------------------------------------------------------------------
	     ! compute expansion
	     !------------------------------------------------------------------
	     !------------------------------------------------------------------
             ! Compute Legendre polynomial box-particle interaction
             !------------------------------------------------------------------
             
             reci      = 1.0_MK/curboxrho(1)
             sine      = SIN(curboxtheta(1))
             cosine    = COS(curboxtheta(1))
             val       = -sine
             prod      = 1.0_MK
             
	     DO m=0,pexp
               Pnm(m,m) = fracfac(m)*prod
               prod     = prod * val
             ENDDO
             
	     DO m=0,pexp-1
                Pnm(m+1,m) = cosine*REAL(2*m + 1,MK)*Pnm(m,m)
             ENDDO
             
	     DO n=2,pexp
                val = cosine*REAL(2*n-1,MK)
                DO m=0,n-1
                   Pnm(n,m)=(val*Pnm(n-1,m)-REAL(n+m-1,MK)* &
                                      Pnm(n-2,m))/REAL(n-m,MK)
                ENDDO
             ENDDO
             
	     
	     !------------------------------------------------------------------
             ! Compute Ynm(n,m) and Ynm(n,-m)
             !------------------------------------------------------------------
	     DO n=0,pexp
                m = 0
                angle = REAL(m,MK)*curboxphi(1)
                Ynm(n,m) = sqrtfac(n,m)*Pnm(n,m)*CMPLX(COS(angle),SIN(angle))
                DO m=1,n
                   angle     = REAL(m,MK)*curboxphi(1)
                   Ynm(n,m)  = sqrtfac(n,m)*Pnm(n,m)*   &
                               CMPLX(COS(angle),SIN(angle))
                   Ynm(n,-m) = CONJG(Ynm(n,m))
                ENDDO
             ENDDO

             !------------------------------------------------------------------
             ! Compute the Outer expansion
             !------------------------------------------------------------------
	     prod = 1.0_MK
             DO n=0,pexp
                prod = prod * curboxrho(1)
                DO m=-n,n
                   Outer(n,m) = (-1)**n*CI**ABS(m)*Ynm(n,m)/(Anm(n,m)*prod)
                ENDDO
             ENDDO

             !------------------------------------------------------------------
             ! Evaluate potential, using multipole expansion coefficients
             !------------------------------------------------------------------
	     DO n=0,pexp
                DO m=-n,n
                   
#if                __DIM == __SFIELD
		   potential(i)=potential(i) + expansion(curbox,n,m) &
                    &           *Outer(n,-m)
#else
                    DO j=1,lda
                       potential(j,i)=potential(j,i) + expansion(j,curbox,n,m) &
                    &           *Outer(n,-m)
                    ENDDO
#endif

                ENDDO
             ENDDO

           ELSE
	   
             !-----------------------------------------------------------------
             !  not far enough, push children if present
             !-----------------------------------------------------------------
	     IF (ppm_debug .GT. 0) THEN
                CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
                &    'not far enough',info)
             ENDIF
             
	     IF (nchld(curbox) .GT. 0) THEN
                !---------------------------------------------------------------
                !  not far enough, push children if present
                !---------------------------------------------------------------
                IF (ppm_debug .GT. 0) THEN
                   CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
                   &    'push children',info)
                ENDIF
               DO j=1,nchld(curbox)
                  stack(stackpointer) = child(j,curbox)
                  stackpointer = stackpointer + 1
               ENDDO
	       
             ELSE
	     
             !------------------------------------------------------------------
             !  no children, direct computation
             !------------------------------------------------------------------
              
	      IF (ppm_debug .GT. 0) THEN
                 CALL ppm_write(ppm_rank,'ppm_fmm_potential', &
                 &    'no children',info)
              ENDIF
              
	      first = lhbx(1,curbox)
              last  = lhbx(2,curbox)
              
	      DO j=first,last !loop over particles in leaf
                  
		 !-------------------------------------------------------------
                 ! Evaluate potential, direct method
                 !-------------------------------------------------------------
                 dx = xp(1,j) - tp(1,i)
                 dy = xp(2,j) - tp(2,i)
                 dz = xp(3,j) - tp(3,i)
                 
                 rad = dx*dx + dy*dy + dz*dz
                 IF(rad.GT.eps)THEN
	            rad = 1.0_MK/SQRT(rad)                  
                    
#if                 __DIM == __SFIELD
		    potential(i) = potential(i) + wp(j)*rad
#else
                    DO l=1,lda
                       potential(l,i) = potential(l,i) + wp(l,j)*rad 
                    ENDDO
#endif
                 ENDIF
              ENDDO
             ENDIF
           ENDIF       
	 ENDDO 
      ENDDO

      
      !-------------------------------------------------------------------------
      !  Nullify data pointers
      !-------------------------------------------------------------------------

      NULLIFY(min_box)
      NULLIFY(max_box)
      NULLIFY(boxcost)
      NULLIFY(centerofbox)
      NULLIFY(radius)

      
      !-------------------------------------------------------------------------
      !  Deallocate local data
      !-------------------------------------------------------------------------

      istat     = 0
      
      ldu1(1)   = 0
      ldu2(1:2) = 0
      CALL ppm_alloc(newlpdx,ldu1,ppm_param_dealloc,info)
      istat = istat + info
      CALL ppm_alloc(stack,ldu1,ppm_param_dealloc,info)
      istat = istat + info
      CALL ppm_alloc(xp,ldu2,ppm_param_dealloc,info)
      istat = istat + info
      CALL ppm_alloc(wp,ldu2,ppm_param_dealloc,info)
      istat = istat + info
      CALL ppm_alloc(part_subtop,ldu2,ppm_param_dealloc,info)
      istat = istat + info
      
      IF (istat .NE. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_fmm_expansion', &
      &       'error deallocating newlpdx',__LINE__,info)
      GOTO 9999
      ENDIF

      
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------

9999  CONTINUE
      
      CALL substop('ppm_fmm_potential',t0,info)
      
      RETURN

#if   (__KIND == __SINGLE_PRECISION && __DIM == __SFIELD)
      END SUBROUTINE ppm_fmm_potential_s_sf
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __SFIELD)
      END SUBROUTINE ppm_fmm_potential_d_sf
#elif (__KIND == __SINGLE_PRECISION && __DIM == __VFIELD)
      END SUBROUTINE ppm_fmm_potential_s_vf
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __VFIELD)
      END SUBROUTINE ppm_fmm_potential_d_vf
#endif

     
  
