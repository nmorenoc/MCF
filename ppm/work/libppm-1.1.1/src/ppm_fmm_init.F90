#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_fmm_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Initialisation of FMM. This routine calls the 
      !                 ppm_tree-routine to get the tree 
      !                 information and stores it. 
      !                 Maps the particles to the leafs of the tree.
      !                 It computes the center of the boxes and the radius
      !                 of the leaf boxes and stores it.
      !
      !  Input        : xp(:,:)      (F) the field points
      !                 wp(:,:)      (F) field particle strenghts
      !                 lda          (I) number of source dimensions
      !                 Nm(:)        (I) number of grid points in the
      !                                  global mesh. (0,0,0) if there is
      !                                  no mesh. If a mesh is present, the
      !                                  box boundaries will be aligned
      !                                  with mesh planes.
      !                 ord          (I) expansion order
      !                 min_dom(:)   (F) the minimum coordinate of the
      !                                  domain
      !                 max_dom(:)   (F) the maximum coordinate of the
      !                                  domain
      !                 maxboxcost   (F) the maximum number of particles
      !                                  allowed in a box
      !  Input/output :     
      !                 Np           (I) the number of field points.
      !
      !  Output       : nrofbox     (I) the total number of all boxes
      !                 info        (I) return status. 0 upon success.
      !
      !  Remarks      :  only useful for freespace boundary conditions
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fmm_init.f,v $
      !  Revision 1.21  2006/09/04 18:34:46  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.20  2006/06/29 10:28:35  pchatela
      !  Added vector strengths support
      !
      !  Revision 1.19  2006/06/20 15:13:35  hiebers
      !  BUGFIX: adjusted indices for ppm_boxid, ppm_subid
      !
      !  Revision 1.18  2006/06/16 07:52:21  hiebers
      !  Added a new list of topo IDs (topoidlist) to prevent overwriting user defined
      !  topologies
      !
      !  Revision 1.17  2005/09/19 13:03:28  polasekb
      !  code cosmetics
      !
      !  Revision 1.16  2005/09/12 13:30:33  polasekb
      !  added ppm_subid
      !
      !  Revision 1.15  2005/09/11 18:05:30  polasekb
      !  (final?) corrected version
      !  (also works parallel :-)
      !
      !  Revision 1.14  2005/09/11 11:43:39  polasekb
      !  moved mapping and second tree call to init
      !
      !  Revision 1.13  2005/08/30 08:48:30  polasekb
      !  added timing for tree
      !
      !  Revision 1.12  2005/08/25 13:51:49  polasekb
      !  corrected data allocation of theta,phi,rho
      !
      !  Revision 1.11  2005/08/11 15:12:53  polasekb
      !  added argument maxboxcost
      !
      !  Revision 1.10  2005/08/08 13:34:25  polasekb
      !  removec fmm_prec
      !
      !  Revision 1.9  2005/08/04 16:00:41  polasekb
      !  moved some allocation to init
      !
      !  Revision 1.8  2005/07/29 12:35:05  polasekb
      !  changed diagonal to radius
      !
      !  Revision 1.7  2005/07/27 14:58:26  polasekb
      !  added new argument wp to subroutine call
      !
      !  Revision 1.6  2005/07/25 15:01:57  polasekb
      !  adapted some tree coefficients
      !
      !  Revision 1.5  2005/07/25 13:39:20  polasekb
      !  bugfix in array indices
      !
      !  Revision 1.4  2005/07/21 13:21:32  polasekb
      !  removed nullify
      !
      !  Revision 1.3  2005/06/02 13:54:55  polasekb
      !  removed totalmass
      !
      !  Revision 1.2  2005/05/30 09:37:24  polasekb
      !  correctet computing of centerofbox
      !
      !  Revision 1.1  2005/05/27 07:53:40  polasekb
      !  Initial Implementation
      !
      !  Revision 0  2004/11/16 15:59:14 polasekb
      !  start
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   (__KIND == __SINGLE_PRECISION && __DIM == __SFIELD)
      SUBROUTINE ppm_fmm_init_s_sf(xp,wp,Np,Nm,ord,min_dom,max_dom,maxboxcost, &
      &          nrofbox,info)
#elif ( __KIND == __DOUBLE_PRECISION && __DIM == __SFIELD)
      SUBROUTINE ppm_fmm_init_d_sf(xp,wp,Np,Nm,ord,min_dom,max_dom,maxboxcost, &
      &          nrofbox,info)
#elif (__KIND == __SINGLE_PRECISION && __DIM == __VFIELD)
      SUBROUTINE ppm_fmm_init_s_vf(xp,wp,lda,Np,Nm,ord,min_dom,max_dom,maxboxcost, &
      &          nrofbox,info)
#elif ( __KIND == __DOUBLE_PRECISION && __DIM == __VFIELD)
      SUBROUTINE ppm_fmm_init_d_vf(xp,wp,lda,Np,Nm,ord,min_dom,max_dom,maxboxcost, &
      &          nrofbox,info)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      
      USE ppm_module_tree
      USE ppm_module_data
      USE ppm_module_data_fmm
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_map
      USE ppm_module_topo_box2subs
      USE ppm_module_topo      
      USE ppm_module_substart
      USE ppm_module_substop 
      USE ppm_module_write
      USE ppm_module_check_topoid
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
#if   __DIM == __SFIELD      
      REAL(MK), DIMENSION(:  ), POINTER       :: wp 
#elif __DIM == __VFIELD
      REAL(MK), DIMENSION(:,:), POINTER       :: wp 
      INTEGER                 , INTENT(IN   ) :: lda
#endif

      REAL(MK), DIMENSION(:,:), POINTER       :: xp
      INTEGER                 , INTENT(INOUT) :: Np
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: Nm
      INTEGER                 , INTENT(IN   ) :: ord
      REAL(MK)                , INTENT(IN   ) :: maxboxcost
      REAL(MK), DIMENSION(:  ), INTENT(IN   ) :: min_dom,max_dom
      INTEGER                 , INTENT(  OUT) :: nrofbox,info

      
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      LOGICAL                             :: pruneboxes,OK
      LOGICAL                             :: TopoExists
      LOGICAL,DIMENSION(3)                :: fixed
      INTEGER,DIMENSION(1)                :: ldu1
      INTEGER,DIMENSION(2)                :: ldu2
      INTEGER,DIMENSION(3)                :: ldu3         
      INTEGER                             :: iopt,i,j,k,l
      CHARACTER(LEN=ppm_char)             :: cbuf
      INTEGER                             :: box,first,last,nrpbox
      REAL(MK)                            :: t0
      REAL(MK)                            :: tmp
      REAL(MK),DIMENSION(3)               :: diagvec
#if   __DIM == __SFIELD
      INTEGER,DIMENSION(3)                :: ldu,ldl
#else 
      INTEGER,DIMENSION(4)                :: ldu,ldl
#endif
      INTEGER                             :: treetype,minboxes
      REAL(MK)                            :: maxvariance
      INTEGER                             :: maxlevels
      REAL(MK),DIMENSION(3,2)             :: weights
      REAL(MK),DIMENSION(3  )             :: minboxsize
      REAL(MK),DIMENSION(:,:), POINTER    :: min_box,max_box
      REAL(MK),DIMENSION(:),   POINTER    :: boxcost       
      REAL(MK),DIMENSION(:)  , POINTER    :: fac,fracfac
      REAL(MK),DIMENSION(:,:), POINTER    :: Anm,sqrtfac
       
      INTEGER                             :: n,m,level
      INTEGER,DIMENSION(:),    POINTER    :: box2proc,boxid 
      REAL(MK),DIMENSION(:),   POINTER    :: cost
      REAL(MK)                            :: ghostsize
      INTEGER                             :: decomp,assig
      INTEGER,DIMENSION(6)                :: bcdef           
      INTEGER,DIMENSION(:),    POINTER    :: subs2proc,isublist
      INTEGER                             :: nsublist 
      REAL(MK),DIMENSION(:,:), POINTER    :: min_sub,max_sub 
      INTEGER                             :: nsubs,topoid  
      INTEGER,DIMENSION(:),    POINTER    :: new_subs2proc
      INTEGER                             :: mapt,Mpart
      INTEGER                             :: istat         
      REAL(MK),DIMENSION(:  ), POINTER    :: radius,totalmass
      REAL(MK),DIMENSION(:,:), POINTER    :: centerofbox      
      REAL(MK),DIMENSION(:,:), POINTER    :: treepart

#if   __DIM == __SFIELD              
      REAL(MK),DIMENSION(:  ), POINTER    :: treewp
#elif __DIM == __VFIELD
      REAL(MK),DIMENSION(:,:), POINTER    :: treewp
#endif

      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------
      CALL substart('ppm_fmm_init',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN  
         DO i=1,ppm_dim
            IF (min_dom(i) .GT. max_dom(i)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_fmm_init',   &
      &               'min_dom must be <= max_dom !',__LINE__,info)
               GOTO 9999
            ENDIF
         ENDDO
      ENDIF

      
      !-------------------------------------------------------------------------
      ! Set init  variables
      !-------------------------------------------------------------------------
      ppm_fmm_initialized        = .TRUE.
      
      !-------------------------------------------------------------------------
      ! Define global variables
      !-------------------------------------------------------------------------
      order        = ord

      !-------------------------------------------------------------------------
      ! Set tree input variables 
      !-------------------------------------------------------------------------
      treetype     = ppm_param_tree_oct
      minboxes     = ppm_nproc
      pruneboxes   = .FALSE.
      minboxsize   = (/0.001_MK,0.001_MK,0.001_MK/)
      maxvariance  = -1.0_MK
      fixed        = (/.FALSE.,.FALSE.,.FALSE./)
      weights(:,1) = (/1.0_MK,0.0_MK,0.0_MK/)
      weights(:,2) = (/0.0_MK,0.0_MK,1.0_MK/)
      maxlevels    = 20
      
      
#if   __KIND == __SINGLE_PRECISION      
      !-------------------------------------------------------------------------
      ! Build the tree - single precision
      !-------------------------------------------------------------------------
      CALL ppm_tree(xp,Np,Nm,min_dom,max_dom,treetype,minboxes,             &
      &        pruneboxes,minboxsize,maxvariance,maxboxcost,maxlevels,fixed,&
      &        weights,min_box_s,max_box_s,lhbx,lpdx,boxcost_s,parent,nchld,&
      &        child,blevel,nbox,nbpl,nlevel,info)
      min_box => min_box_s
      max_box => max_box_s
      boxcost => boxcost_s

      
#else
      !-------------------------------------------------------------------------
      ! Build the tree - double precision
      !-------------------------------------------------------------------------
      CALL ppm_tree(xp,Np,Nm,min_dom,max_dom,treetype,minboxes,             &
      &        pruneboxes,minboxsize,maxvariance,maxboxcost,maxlevels,fixed,&
      &        weights,min_box_d,max_box_d,lhbx,lpdx,boxcost_d,parent,nchld,&
      &        child,blevel,nbox,nbpl,nlevel,info)
      min_box => min_box_d
      max_box => max_box_d
      boxcost => boxcost_d
#endif

      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_init',  &
       &       'Calling tree failed.',__LINE__,info)
          GOTO 9999
      ENDIF
      
      nrofbox = nbox
      
      IF (ppm_debug.GT.0) THEN  
         CALL ppm_write(ppm_rank,'ppm_fmm_init', &
         &    'calling tree successful',info)
         WRITE (cbuf,'(A,I8)') 'nbox = ',nbox
         CALL ppm_write(ppm_rank,'ppm_fmm_init',cbuf,info)
      ENDIF

      !-------------------------------------------------------------------------
      ! Allocate topoidlist
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu1(1) = nlevel
      CALL ppm_alloc(topoidlist,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating topoidlist',__LINE__,info)
      GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      ! Allocate boxpart
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu1(1) = Np
      CALL ppm_alloc(boxpart,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating boxpart',__LINE__,info)
      GOTO 9999
      ENDIF
      DO i=1,Np
         boxpart(i) = 0.0_MK
      ENDDO
      
      !-------------------------------------------------------------------------
      ! Store which particle is in which leaf box
      !-------------------------------------------------------------------------
      IF (Np .GT. 0) THEN
        DO i=1,nbox
          IF (nchld(i) .EQ. 0) THEN
             DO j=lhbx(1,i),lhbx(2,i)
                boxpart(lpdx(j)) = i
             ENDDO
           ENDIF
        ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      ! make top level topology seperate
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
      ! Allocate ppm_boxid, ppm_subid, box2proc and cost
      !-------------------------------------------------------------------------

      ldu2(1) = nbox
      ldu2(2) = nlevel
      CALL ppm_alloc(ppm_boxid,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating ppm_boxid',__LINE__,info)
         GOTO 9999
      ENDIF
      
      CALL ppm_alloc(ppm_subid,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating ppm_boxid',__LINE__,info)
         GOTO 9999
      ENDIF

      ldu1(1) = nbox
      CALL ppm_alloc(box2proc,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating box2proc',__LINE__,info)
         GOTO 9999
      ENDIF

      ldu1(1) = MAXVAL(nbpl(:))
      CALL ppm_alloc(cost,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating box2proc',__LINE__,info)
         GOTO 9999
      ENDIF
      
      
      !-------------------------------------------------------------------------
      ! Set parallelisation input variables
      !-------------------------------------------------------------------------
      DO i=1,nbox
         box2proc(i)       = 0.0_MK
         DO j=1,nlevel
            ppm_boxid(i,j) = 0.0_MK
            ppm_subid(i,j) = 0.0_MK
         ENDDO
      ENDDO

      DO i=1,size(cost,1)
         cost(i)      = 1.0_MK
      ENDDO

      decomp     = ppm_param_decomp_user_defined
      assig      = ppm_param_assign_internal
      ghostsize  = 0.0_MK
      bcdef(1:6) = ppm_param_bcdef_freespace
      
      
      !-------------------------------------------------------------------------
      ! Transforming leaf boxes into subdomains and create the topologies
      !-------------------------------------------------------------------------
      
      CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,max_sub, &
           &                 nsubs,info,boxid,-level,blevel,child)
      IF (info.NE.0) THEN
         CALL ppm_write(ppm_rank,'ppm_fmm_init', &
         &    'topo_box2subs failed',info)
      ENDIF


      
      !-------------------------------------------------------------------------
      ! Find first topo id that is available 
      !-------------------------------------------------------------------------
      k = 1
      CALL ppm_check_topoid(ppm_param_id_user,k,TopoExists,info)
      DO WHILE(TopoExists)
         k = k+1
         CALL ppm_check_topoid(ppm_param_id_user,k,TopoExists,info)
      ENDDO

      topoid = k
      topoidlist(level) = k

      
      !-------------------------------------------------------------------------
      ! Create first topology based on leaf boxes
      !-------------------------------------------------------------------------
      CALL ppm_mktopo(decomp,assig,min_dom,max_dom,bcdef,ghostsize, &
           &          topoid,min_sub,max_sub,cost,subs2proc,nsubs,isublist, &
           &          nsublist,info)
      IF (info.NE.0) THEN
         CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_init', &
              &         'mktopo failed',__LINE__,info)
      ENDIF

      DO i = 1,nsubs
         ppm_boxid(i,level) = boxid(i)
      ENDDO
      DO j=1,SIZE(boxid)
         ppm_subid(boxid(j),level) = j
      ENDDO
      DO i=1,nsubs
         box2proc(boxid(i)) = subs2proc(i)
      ENDDO

      
      !-------------------------------------------------------------------------
      ! Loop over the levels of the tree and register each level as topology
      !-------------------------------------------------------------------------
      assig      = ppm_param_assign_user_defined
      
      DO i=level+1,nlevel
        !-----------------------------------------------------------------------
        ! Assigning topoids  that are not used before
        !-----------------------------------------------------------------------
         
         CALL ppm_check_topoid(ppm_param_id_user,k,TopoExists,info)
         DO WHILE(TopoExists)
            k = k+1
            CALL ppm_check_topoid(ppm_param_id_user,k,TopoExists,info)
         ENDDO
         topoid = k
         topoidlist(i) = k
         k = k+1


        !-----------------------------------------------------------------------
        ! Call subroutine to get subs
        !-----------------------------------------------------------------------

        CALL ppm_topo_box2subs(min_box,max_box,nchld,nbox,min_sub,max_sub, &
          &                      nsubs,info,boxid,-topoid,blevel,child)
        IF (info.NE.0) THEN
           CALL ppm_write(ppm_rank,'ppm_fmm_init', &
           &    'topo_box2subs failed',info)
        ENDIF

	
        !-----------------------------------------------------------------------
        ! Allocate new subs2proc
        !-----------------------------------------------------------------------
        iopt = ppm_param_alloc_grow
        ldu1(1) = nsubs
        CALL ppm_alloc(new_subs2proc,ldu1,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
        &       'error allocating new_subs2proc',__LINE__,info)
        GOTO 9999
        ENDIF

        DO j=1,nsubs
           new_subs2proc(j) = box2proc(parent(boxid(j)))
           box2proc(boxid(j)) = new_subs2proc(j)
        ENDDO

	
        !-----------------------------------------------------------------------
        ! Call ppm_mktopo to get topology
        !-----------------------------------------------------------------------
	CALL ppm_mktopo(decomp,assig,min_dom,max_dom,bcdef,ghostsize, & 
        &          topoid,min_sub,max_sub,cost,new_subs2proc,nsubs, &  
        &          isublist,nsublist,info)
        IF (info.NE.0) THEN
           CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_init', &
           &         'mktopo failed',__LINE__,info)
        ENDIF
        ppm_boxid(1:nsubs,i) = boxid(1:nsubs)
        DO j=1,SIZE(boxid)
           ppm_subid(boxid(j),i)   = j
        ENDDO
      ENDDO 

      
      !-------------------------------------------------------------------------
      ! Map for the lowest level (leafs of tree)
      !-------------------------------------------------------------------------
      topoid = topoidlist(nlevel)

      IF (ppm_nproc .GT. 1) THEN

      !-------------------------------------------------------------------------
      !  Map the particles onto the finest tree topology = topoid=nlevel
      !-------------------------------------------------------------------------
      mapt = ppm_param_map_global
      CALL ppm_map_part(xp,ppm_dim,Np,Mpart,topoid,mapt,info)   ! positions
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to start global mapping.',info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Push along the strength of the particles and the boxpart
      !-------------------------------------------------------------------------      
      mapt = ppm_param_map_push
#if   __DIM == __SFIELD      
      CALL ppm_map_part(wp,Np,Mpart,topoid,mapt,info)   ! strengths
#else
      CALL ppm_map_part(wp,lda,Np,Mpart,topoid,mapt,info)   ! strengths      
#endif
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF
      
      CALL ppm_map_part(boxpart,Np,Mpart,topoid,mapt,info)   ! boxpart
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF
      
      mapt = ppm_param_map_send
#if   __DIM == __SFIELD      
      CALL ppm_map_part(wp,Np,Mpart,topoid,mapt,info)   ! strengths
#else
      CALL ppm_map_part(wp,lda,Np,Mpart,topoid,mapt,info)   ! strengths      
#endif
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to send particles.',info)
          GOTO 9999
      ENDIF
      
      mapt = ppm_param_map_pop
      CALL ppm_map_part(boxpart,Np,Mpart,topoid,mapt,info)   ! boxpart
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to push strengths.',info)
          GOTO 9999
      ENDIF
      
#if   __DIM == __SFIELD      
      CALL ppm_map_part(wp,Np,Mpart,topoid,mapt,info)   ! strengths
#else
      CALL ppm_map_part(wp,lda,Np,Mpart,topoid,mapt,info)   ! strengths      
#endif
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to pop strengths.',info)
          GOTO 9999
      ENDIF
      
      CALL ppm_map_part(xp,ppm_dim,Np,Mpart,topoid,mapt,info)   ! positions
      IF (info .NE. 0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Failed to pop positions.',info)
          GOTO 9999
      ENDIF
                                                                     
      !-------------------------------------------------------------------------
      !  Update and store and new number of particles 
      !-------------------------------------------------------------------------
      
      Np = Mpart
      
      IF (ppm_debug .GT. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_fmm_init', &
        &             'Done mapping particles.',info)
        WRITE(cbuf,'(A,I8)') 'Local number of particles is now: ',Np
        CALL ppm_write(ppm_rank,'ppm_fmm_init',cbuf,info)
      ENDIF 
      
                                                                     
      !-------------------------------------------------------------------------
      !  Check that particles have been mapped correctly
      !-------------------------------------------------------------------------
      CALL ppm_topo_check(xp,Np,OK,info)
      IF (info .NE. 0) THEN
         CALL ppm_write(ppm_rank,'ppm_fmm_init', &
         &    'Failed to check topology.',info)
      ENDIF
      
      IF (.NOT.OK) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'Particles not mapped correctly!',info) 
          GOTO 9999
      ENDIF
      
      
      !-------------------------------------------------------------------------
      ! Rebuild the tree to get the correct lpdx and lhbx arrays
      ! (after mapping particles changed)
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      !-------------------------------------------------------------------------
      ! Build the tree - single precision
      !-------------------------------------------------------------------------
      CALL ppm_tree(xp,Np,Nm,min_dom,max_dom,treetype,minboxes,        &
      &        pruneboxes,minboxsize,maxvariance,maxboxcost,maxlevels,fixed,&
      &        weights,min_box_s,max_box_s,lhbx,lpdx,boxcost_s,parent,nchld,&
      &        child,blevel,nbox,nbpl,nlevel,info)
      min_box => min_box_s
      max_box => max_box_s
      boxcost => boxcost_s

      
#else
      !-------------------------------------------------------------------------
      ! Build the tree - double precision
      !-------------------------------------------------------------------------
      CALL ppm_tree(xp,Np,Nm,min_dom,max_dom,treetype,minboxes,        &
      &        pruneboxes,minboxsize,maxvariance,maxboxcost,maxlevels,fixed,&
      &        weights,min_box_d,max_box_d,lhbx,lpdx,boxcost_d,parent,nchld,&
      &        child,blevel,nbox,nbpl,nlevel,info)
      min_box => min_box_d
      max_box => max_box_d
      boxcost => boxcost_d
#endif

      IF (info .NE. 0) THEN
          info = ppm_error_error
          CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_init',  &
          &    'Calling tree (2) failed.',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      ! Store number of boxes in nrofbox
      !-------------------------------------------------------------------------
      nrofbox = nbox

      IF (ppm_debug.GT.0) THEN
          CALL ppm_write(ppm_rank,'ppm_fmm_init', &
          &    'calling tree (2) successful',info)
          WRITE (cbuf,'(A,I8)') 'nbox = ',nbox
          CALL ppm_write(ppm_rank,'ppm_fmm_init',cbuf,info)
      ENDIF

      ENDIF ! ppm_nrpco .GT. 1

       
      !-------------------------------------------------------------------------
      ! Allocate all data (single/double)
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! Allocate and initialise sqrtfac,Anm, Outer, fracfac (single/double)
      !-------------------------------------------------------------------------

      iopt = ppm_param_alloc_fit

#if   __KIND == __SINGLE_PRECISION
      ldl(1) = 0
      ldl(2) = -order
      ldu(1) = order
      ldu(2) = order
      CALL ppm_alloc(sqrtfac_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating sqrtfac',__LINE__,info)
         GOTO 9999
      ENDIF
      
      CALL ppm_alloc(Anm_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Anm',__LINE__,info)
         GOTO 9999
      ENDIF
      
      CALL ppm_alloc(Outer_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Outer',__LINE__,info)
         GOTO 9999
      ENDIF
      
      ldl(1) = 0
      ldu(1) = order
      CALL ppm_alloc(fracfac_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating frac_fac',__LINE__,info)
         GOTO 9999
      ENDIF
      
      ! Initiatise variables
      Anm     => Anm_s
      sqrtfac => sqrtfac_s
      fracfac => fracfac_s


      
#else
      ldl(1) = 0
      ldl(2) = -order
      ldu(1) = order
      ldu(2) = order
      CALL ppm_alloc(sqrtfac_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating sqrtfac',__LINE__,info)
      GOTO 9999
      ENDIF
      
      CALL ppm_alloc(Anm_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Anm',__LINE__,info)
      GOTO 9999
      ENDIF
      
      CALL ppm_alloc(Outer_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Outer',__LINE__,info)
      GOTO 9999
      ENDIF
      
      ldl(1) = 0
      ldu(1) = order
      CALL ppm_alloc(fracfac_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating fracfac',__LINE__,info)
      GOTO 9999
      ENDIF

      
      ! Initiatise variables
      Anm     => Anm_d
      sqrtfac => sqrtfac_d
      fracfac => fracfac_d

#endif








      !-------------------------------------------------------------------------
      ! Allocate and initialise expansions
      !-------------------------------------------------------------------------
#if   (__KIND == __SINGLE_PRECISION && __DIM == __SFIELD)
      ldl(1)  = 1
      ldl(2)  = 0
      ldl(3)  = -order
      ldu(1)  = nbox
      ldu(2)  = order
      ldu(3)  = order
      CALL ppm_alloc(expansion_s_sf,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating expansion',__LINE__,info)
         GOTO 9999
      ENDIF  
      
      DO j=0,order
         DO k=-order,order
            DO i=1,nbox
               expansion_s_sf(i,j,k) = 0.0_MK
            ENDDO
         ENDDO
      ENDDO
 
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __SFIELD)      
      ldl(1)  = 1
      ldl(2)  = 0
      ldl(3)  = -order
      ldu(1)  = nbox
      ldu(2)  = order
      ldu(3)  = order
      CALL ppm_alloc(expansion_d_sf,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating expansion',__LINE__,info)
         GOTO 9999
      ENDIF

      DO j=0,order
         DO k=-order,order
            DO i=1,nbox
               expansion_d_sf(i,j,k) = 0.0_MK
            ENDDO
         ENDDO
      ENDDO
      
#elif (__KIND == __SINGLE_PRECISION && __DIM == __VFIELD)
      ldl(1)  = 1
      ldl(2)  = 1
      ldl(3)  = 0
      ldl(4)  = -order
      ldu(1)  = lda
      ldu(2)  = nbox
      ldu(3)  = order
      ldu(4)  = order
      CALL ppm_alloc(expansion_s_vf,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating expansion',__LINE__,info)
         GOTO 9999
      ENDIF

      DO j=0,order
         DO k=-order,order
            DO i=1,nbox
               DO l=1,lda
                  expansion_s_vf(l,i,j,k) = 0.0_MK
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __VFIELD)      
      ldl(1)  = 1
      ldl(2)  = 1
      ldl(3)  = 0
      ldl(4)  = -order
      ldu(1)  = lda
      ldu(2)  = nbox
      ldu(3)  = order
      ldu(4)  = order   
      CALL ppm_alloc(expansion_d_vf,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating expansion',__LINE__,info)
         GOTO 9999
      ENDIF

      DO j=0,order
         DO k=-order,order
            DO i=1,nbox
               DO l=1,lda
                  expansion_d_vf(l,i,j,k) = 0.0_MK
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
#endif      
      

      !-------------------------------------------------------------------------
      ! Allocate multipole coefficient variables      
      ! Pnm, Ynm, fac, rho, theta, phi      
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      ldl(1) = 0
      ldu(1) = 2*order
      CALL ppm_alloc(fac_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating fac',__LINE__,info)
         GOTO 9999
      ENDIF
      
      ldu1(1) = ppm_nproc*Np
      CALL ppm_alloc(rho_s,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating rho',__LINE__,info)
         GOTO 9999
      ENDIF
      
      CALL ppm_alloc(theta_s,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating theta',__LINE__,info)
         GOTO 9999
      ENDIF
      
      CALL ppm_alloc(phi_s,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating phi',__LINE__,info)
         GOTO 9999
      ENDIF

      DO i = 1,ldu1(1)
         rho_s(i)   = 0.0_MK
         theta_s(i) = 0.0_MK
         phi_s(i)   = 0.0_MK
      ENDDO
      
      ldl(1) = 0
      ldl(2) = -order
      ldu(1) = order
      ldu(2) = order
      CALL ppm_alloc(Ynm_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating Ynm',__LINE__,info)
      GOTO 9999
      ENDIF
      
      CALL ppm_alloc(Pnm_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating Pnm',__LINE__,info)
         GOTO 9999
      ENDIF
      
      ldl(1) = 0
      ldl(2) = -2*order
      ldu(1) = 2*order
      ldu(2) = 2*order
      CALL ppm_alloc(Inner_s,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating Inner',__LINE__,info)
         GOTO 9999
      ENDIF
      
      DO i=ldl(1),ldu(1)
         DO j=ldl(2),ldu(2)
            Inner_s(i,j) = 0.0_MK
         ENDDO
      ENDDO
      ! Initiatise variables
      fac     => fac_s
#else
      ldl(1) = 0
      ldu(1) = 2*order
      CALL ppm_alloc(fac_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating fac',__LINE__,info)
         GOTO 9999
      ENDIF
      
      ldu1(1) = ppm_nproc*Np
      CALL ppm_alloc(rho_d,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating rho',__LINE__,info)
        GOTO 9999
      ENDIF
      
      CALL ppm_alloc(theta_d,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating theta',__LINE__,info)
         GOTO 9999
      ENDIF
      
      CALL ppm_alloc(phi_d,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating phi',__LINE__,info)
         GOTO 9999
      ENDIF

      DO i = 1,ldu1(1)
         rho_d(i)   = 0.0_MK
         theta_d(i) = 0.0_MK
         phi_d(i)   = 0.0_MK
      ENDDO
      
      ldl(1) = 0
      ldl(2) = -order
      ldu(1) = order
      ldu(2) = order
      CALL ppm_alloc(Ynm_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating Ynm',__LINE__,info)
         GOTO 9999
      ENDIF
      
      CALL ppm_alloc(Pnm_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating Pnm',__LINE__,info)
         GOTO 9999
      ENDIF
      
      ldl(1) = 0
      ldl(2) = -2*order
      ldu(1) = 2*order
      ldu(2) = 2*order
      CALL ppm_alloc(Inner_d,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating Inner',__LINE__,info)
         GOTO 9999
      ENDIF
      
      
      DO i=ldl(1),ldu(1)
         DO j=ldl(2),ldu(2)
            Inner_d(i,j) = 0.0_MK
         ENDDO
      ENDDO

      
      ! Initiatise variables
      fac     => fac_d


#endif
 

      !-------------------------------------------------------------------------
      ! Initialise fac with zero
      !-------------------------------------------------------------------------
      DO i = 0,2*order
         fac(i)   = 0.0_MK
      ENDDO
   
      !-------------------------------------------------------------------------
      ! Compute fac, rho, phi, theta and topoid
      !-------------------------------------------------------------------------
      fac(0) = 1
      
      DO i=1,order*2
         fac(i) = fac(i-1)*REAL(i,MK)
      ENDDO
      
      DO n=0,order
         DO m=-n,n
            sqrtfac(n,m) = SQRT(fac(n-ABS(m))/fac(n+ABS(m)))
            Anm(n,m)     = (-1)**n/SQRT(fac(n-m)*fac(n+m))
         ENDDO
      ENDDO

      DO m=0,order
         fracfac(m) = fac(2*m)/(2.**m*fac(m))
      ENDDO
      





      !-------------------------------------------------------------------------
      ! Allocate Cnm
      !-------------------------------------------------------------------------
#if   (__KIND == __SINGLE_PRECISION && __DIM == __SFIELD)      
      ldl(1) = 0
      ldl(2) = -order
      ldu(1) = order
      ldu(2) = order
      CALL ppm_alloc(Cnm_s_sf,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating Cnm',__LINE__,info)
         GOTO 9999
      ENDIF

#elif (__KIND == __DOUBLE_PRECISION && __DIM == __SFIELD)      
      ldl(1) = 0
      ldl(2) = -order
      ldu(1) = order
      ldu(2) = order
      CALL ppm_alloc(Cnm_d_sf,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating Cnm',__LINE__,info)
         GOTO 9999
      ENDIF
      
#elif (__KIND == SINGLE_PRECISION && __DIM == __VFIELD) 
      ldl(1) = 1
      ldl(2) = 0
      ldl(3) = -order
      ldu(1) = lda
      ldu(2) = order
      ldu(3) = order
      CALL ppm_alloc(Cnm_s_vf,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating Cnm',__LINE__,info)
         GOTO 9999
      ENDIF
       
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __VFIELD)      
      ldl(1) = 1
      ldl(2) = 0
      ldl(3) = -order
      ldu(1) = lda
      ldu(2) = order
      ldu(3) = order
      CALL ppm_alloc(Cnm_d_vf,ldl,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating Cnm',__LINE__,info)
         GOTO 9999
      ENDIF

#endif      

    
      !-------------------------------------------------------------------------
      ! Allocate radius (single/double)
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit

#if   __KIND == __SINGLE_PRECISION
      ldu1 = nbox
      CALL ppm_alloc(radius_s,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating radius',__LINE__,info)
         GOTO 9999
      ENDIF
      
#else
      ldu1 = nbox
      CALL ppm_alloc(radius_d,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating radius',__LINE__,info)
         GOTO 9999
      ENDIF

#endif


      !-------------------------------------------------------------------------
      ! Allocate centerofbox (single/double)
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      ldu2(1) = ppm_dim
      ldu2(2) = nbox
      CALL ppm_alloc(centerofbox_s,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating centerofbox',__LINE__,info)
         GOTO 9999
      ENDIF

#else
      ldu2(1) = ppm_dim
      ldu2(2) = nbox
      CALL ppm_alloc(centerofbox_d,ldu2,iopt,info)
      IF (info .NE. 0) THEN 
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
          &    'error allocating centerofbox',__LINE__,info)
          GOTO 9999
      ENDIF
      
#endif


      !-------------------------------------------------------------------------
      ! Allocate totalmass (single/double)
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      ldu1 = nbox
      CALL ppm_alloc(totalmass_s,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
         &    'error allocating totalmass',__LINE__,info)
         GOTO 9999
      ENDIF
      
#else
      ldu1 = nbox
      CALL ppm_alloc(totalmass_d,ldu1,iopt,info)
      IF (info .NE. 0) THEN 
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
          &    'error allocating totalmass',__LINE__,info)
          GOTO 9999
      ENDIF

#endif

      IF (ppm_debug.GT.0) THEN  
         CALL ppm_write(ppm_rank,'ppm_fmm_init','alloc data successful' &
         &    ,info)
      ENDIF


      !-------------------------------------------------------------------------
      ! Check precision and pointing to the correct variables
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      centerofbox  => centerofbox_s
      totalmass    => totalmass_s
      radius       => radius_s
      maxboxcost_s = maxboxcost
#else
      centerofbox  => centerofbox_d
      totalmass    => totalmass_d
      radius       => radius_d
      maxboxcost_d = maxboxcost
#endif

      ! initialise centerofbox with zero

      DO j=1,nbox
         DO i=1,ppm_dim
            centerofbox(i,j) = 0.0_MK
         ENDDO
      ENDDO


      !-------------------------------------------------------------------------
      ! Allocate treepart and treewp
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit

      ldu2(1) = ppm_dim
      ldu2(2) = Np
      CALL ppm_alloc(treepart,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating treepart',__LINE__,info)
         GOTO 9999
      ENDIF

#if   __DIM == __SFIELD            
      ldu1 = Np
      CALL ppm_alloc(treewp,ldu1,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating treewp',__LINE__,info)
         GOTO 9999
      ENDIF

#elif __DIM == __VFIELD
      ldu2(1) = lda
      ldu2(2) = Np
      CALL ppm_alloc(treewp,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_init', &
      &       'error allocating treewp',__LINE__,info)
         GOTO 9999
      ENDIF  

#endif      

      !-------------------------------------------------------------------------
      ! Loop over the boxes to compute the centerofbox of each
      !-------------------------------------------------------------------------
      DO box=1,nbox
         IF (nchld(box) .NE. 0) CYCLE
         first = lhbx(1,box)
         last  = lhbx(2,box)
	 
	 
         !----------------------------------------------------------------------
         ! Compute new array with particle order from tree
	 ! and necessary information (nr. of part. in box)
         !----------------------------------------------------------------------
         nrpbox = last-first+1
         
         
	 DO j=first,last
            treepart(1,j) = xp(1,lpdx(j))
            treepart(2,j) = xp(2,lpdx(j))
            IF(ppm_dim.EQ.3)THEN
               treepart(3,j) = xp(3,lpdx(j))
            ENDIF

#if         __DIM == __SFIELD        
            treewp(j)   = wp(lpdx(j))

#else
            DO i=1,lda
               treewp(i,j) = wp(i,lpdx(j))
            ENDDO
#endif 
         ENDDO


         !----------------------------------------------------------------------
         ! Computing the total mass
         !----------------------------------------------------------------------

#if      __DIM == __SFIELD         
         totalmass(box) = 0.0_MK
         DO i=first,last
            totalmass(box) = totalmass(box) + ABS(treewp(i))
         ENDDO
#else 
         totalmass(box) = 0.0_MK
         DO i=first,last
            DO j=1,lda
               totalmass(box) = totalmass(box) + treewp(j,i)*treewp(j,i)
            ENDDO
         ENDDO 
         totalmass(box) = SQRT(totalmass(box))
#endif 


         !----------------------------------------------------------------------
         ! Compute the centers of the leaf boxes
         !----------------------------------------------------------------------
         
         IF (nrpbox .GT. 0) THEN

#if         __DIM == __SFIELD
            IF(ppm_dim.EQ.2)THEN
               centerofbox(1,box) = 0.0_MK
               centerofbox(2,box) = 0.0_MK
               DO j=first,last
                  centerofbox(1,box) = centerofbox(1,box) + treepart(1,j)* &
                              &   ABS(treewp(j))
                  centerofbox(2,box) = centerofbox(2,box) + treepart(2,j)* &
                              &   ABS(treewp(j))
               ENDDO
               tmp = 1.0_MK/totalmass(box)
               centerofbox(1,box) = centerofbox(1,box) * tmp
               centerofbox(2,box) = centerofbox(2,box) * tmp
            ENDIF
            IF(ppm_dim.EQ.3)THEN
               centerofbox(1,box) = 0.0_MK
               centerofbox(2,box) = 0.0_MK
               centerofbox(3,box) = 0.0_MK
               DO j=first,last
                  centerofbox(1,box) = centerofbox(1,box) + treepart(1,j)* &
                              &   ABS(treewp(j))
                  centerofbox(2,box) = centerofbox(2,box) + treepart(2,j)* &
                              &   ABS(treewp(j))
                  centerofbox(3,box) = centerofbox(3,box) + treepart(3,j)* &
                              &   ABS(treewp(j))
               ENDDO
               tmp = 1.0_MK/totalmass(box)
               centerofbox(1,box) = centerofbox(1,box) *tmp
               centerofbox(2,box) = centerofbox(2,box) *tmp
               centerofbox(3,box) = centerofbox(3,box) *tmp
            ENDIF

#else
            IF(ppm_dim.EQ.2)THEN

               centerofbox(1,box) = 0.0_MK
               centerofbox(2,box) = 0.0_MK
               DO j=first,last
                  tmp = 0.0_MK
                  DO l = 1,lda
                    tmp = tmp + treewp(l,j)*treewp(l,j)
                  ENDDO
                  centerofbox(1,box) = centerofbox(1,box)+ &
                  &          SQRT(tmp)*treepart(1,j)
                  centerofbox(2,box) = centerofbox(2,box)+ &
                  &          SQRT(tmp)*treepart(2,j)
               ENDDO
               tmp = 1.0_MK/totalmass(box)
               centerofbox(1,box) = centerofbox(1,box) * tmp
               centerofbox(2,box) = centerofbox(2,box) * tmp

            ENDIF
            IF(ppm_dim.EQ.3)THEN

               centerofbox(1,box) = 0.0_MK
               centerofbox(2,box) = 0.0_MK
               centerofbox(3,box) = 0.0_MK
               DO j=first,last
                  tmp = 0.0_MK
                  DO l = 1,lda
                    tmp = tmp + treewp(l,j)*treewp(l,j)
                  ENDDO
                  centerofbox(1,box) = centerofbox(1,box)+ &
                  &          SQRT(tmp)*treepart(1,j)
                  centerofbox(2,box) = centerofbox(2,box)+ &
                  &          SQRT(tmp)*treepart(2,j)
                  centerofbox(3,box) = centerofbox(3,box)+ &
                  &          SQRT(tmp)*treepart(3,j)
               ENDDO
               tmp = 1.0_MK/totalmass(box)
               centerofbox(1,box) = centerofbox(1,box) * tmp
               centerofbox(2,box) = centerofbox(2,box) * tmp
               centerofbox(3,box) = centerofbox(3,box) * tmp

            ENDIF


#endif
         
	 ELSE

            IF(ppm_dim.EQ.2)THEN
               centerofbox(1,box) = 0.5_MK*(max_box(1,box) + min_box(1,box))
               centerofbox(2,box) = 0.5_MK*(max_box(2,box) + min_box(2,box))
            ENDIF
            IF(ppm_dim.EQ.3)THEN
               centerofbox(1,box) = 0.5_MK*(max_box(1,box) + min_box(1,box))
               centerofbox(2,box) = 0.5_MK*(max_box(2,box) + min_box(2,box))
               centerofbox(3,box) = 0.5_MK*(max_box(3,box) + min_box(3,box))
            ENDIF

         ENDIF
      ENDDO


      IF (ppm_debug.GT.0) THEN  
         CALL ppm_write(ppm_rank,'ppm_fmm_init','computed centers',info)
      ENDIF

      
      !-------------------------------------------------------------------------
      ! Compute the radius of the leaf boxes
      !-------------------------------------------------------------------------
      DO box=1,nbox
         radius(box) = -1.0_MK
      ENDDO

      IF(ppm_dim.EQ.2)THEN
         DO box=1,nbox
            IF (nchld(box) .NE. 0) CYCLE
            first = lhbx(1,box)
            last  = lhbx(2,box)

            DO j=first,last
               diagvec(1) = treepart(1,j) - centerofbox(1,box)
               diagvec(2) = treepart(2,j) - centerofbox(2,box)
               tmp = diagvec(1)*diagvec(1) + diagvec(2)*diagvec(2)
               tmp = SQRT(tmp)
               IF (tmp .GT. radius(box)) THEN
                  radius(box) = tmp
               ENDIF  
            ENDDO
         ENDDO

      ENDIF
      IF(ppm_dim.EQ.3)THEN
         DO box=1,nbox
            IF (nchld(box) .NE. 0) CYCLE
            first = lhbx(1,box)
            last  = lhbx(2,box)

            DO j=first,last
               diagvec(1) = treepart(1,j) - centerofbox(1,box)
               diagvec(2) = treepart(2,j) - centerofbox(2,box)
               diagvec(3) = treepart(3,j) - centerofbox(3,box)
               tmp = diagvec(1)*diagvec(1) + diagvec(2)*diagvec(2)
               tmp = tmp                   + diagvec(3)*diagvec(3)
               tmp = SQRT(tmp)
               IF (tmp .GT. radius(box)) THEN
                  radius(box) = tmp
               ENDIF  
            ENDDO
         ENDDO

      ENDIF

      IF (ppm_debug.GT.0) THEN  
         CALL ppm_write(ppm_rank,'ppm_fmm_init','computed radius',info)
      ENDIF

      !-------------------------------------------------------------------------
      !  deallocate local variables
      !-------------------------------------------------------------------------
      ldu1    = 0
      ldu2(1) = 0
      ldu2(2) = 0
      istat   = 0
      iopt = ppm_param_dealloc
      CALL ppm_alloc(treepart,ldu2,iopt,info)
      istat=istat+info
      CALL ppm_alloc(treewp,ldu1,iopt,info)
      istat=istat+info
      CALL ppm_alloc(box2proc,ldu1,iopt,info)
      istat=istat+info
      CALL ppm_alloc(subs2proc,ldu1,iopt,info)
      istat=istat+info
      CALL ppm_alloc(isublist,ldu1,iopt,info)
      istat=istat+info
      CALL ppm_alloc(boxid,ldu1,iopt,info)
      istat=istat+info
      CALL ppm_alloc(cost,ldu1,iopt,info)
      istat=istat+info
      CALL ppm_alloc(new_subs2proc,ldu1,iopt,info)
      istat=istat+info
      CALL ppm_alloc(subs2proc,ldu1,iopt,info)
      istat=istat+info

      IF (istat .NE. 0) THEN
          WRITE(cbuf,'(A,I3,A)') 'for ',istat,'error while dealloc'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fmm_init',cbuf,__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Nullify data pointers
      !-------------------------------------------------------------------------
      NULLIFY(min_box)
      NULLIFY(max_box)
      NULLIFY(boxcost)
      NULLIFY(centerofbox)
      NULLIFY(totalmass)
      NULLIFY(radius)

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
9999  CONTINUE
      CALL substop('ppm_fmm_init',t0,info)
      RETURN

#if   (__KIND == __SINGLE_PRECISION && __DIM == __SFIELD)
      END SUBROUTINE ppm_fmm_init_s_sf
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __SFIELD)
      END SUBROUTINE ppm_fmm_init_d_sf
#elif (__KIND == __SINGLE_PRECISION && __DIM == __VFIELD)
      END SUBROUTINE ppm_fmm_init_s_vf
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __VFIELD)
      END SUBROUTINE ppm_fmm_init_d_vf
#endif
      
