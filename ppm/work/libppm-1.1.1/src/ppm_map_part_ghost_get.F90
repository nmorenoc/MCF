        !------------------------------------------------------------------------
        ! Subroutine  :  ppm_map_part_ghost_get
        !------------------------------------------------------------------------
        !
        ! Purpose     : This routine maps/adds the ghost particles on the 
        !               current topology. This routine is similar to the 
        !               partial mapping routine (ppm_map_part_partial) in the 
        !               sense that the ghost particles are assumed to be 
        !               located on neighbouring processors only, and thus only
        !               require a nearest neighbour communication. 
        ! 
        ! Input       :  xp(:,:)      (F) : the position of the particles
        !                lda          (I) : leading dimension of xp
        !                Npart        (I) : the number of particles (on the
        !                                    local processor)
        !                isymm        (I) : indicator for the use of symmetry 
        !                                    isymm > 0 use symmetry
        !                                    isymm = 0 do not use symmetry
        !                ghostsize    (F) : the size of the ghost layer
        !         
        !                shear_length (F) : the length sheard in each side of
        !                                    the images for Lees-Edwards
        !                                    boundaries.
        !                                    
        !                                     shear_length(direct,dim)
        !                                     dim   :same index as in bcdef(dim)
        !                                     direct: 1, x ; 2, y; 3, z.
        !                                     2D : dim from 1 to 4
        !                                          direct from 1 to 2 
        !                                     3D : dim from 1 to 6
        !                                          direct from 1 to 3 
        !                                   
        !                                    
        ! Input/output: info         (I) : return status, 0 on success
        !
        ! Remarks     :  This routine SHOULD be efficient - since it will be 
        !                called frequently. One way of improving (?) the
        !                performance is to use the cell lists to find the 
        !                potential ghosts.
        !
        !                This routine should be split in several routines !
        !
        !                The ghosts are found in three steps: 
        !
        !                   1) consider ghosts within the local processor
        !                   2) consider ghosts on neighbouring processors
        !                   1) consider ghosts on periodic images of neighbouring
        !                      processors
        !                  
        !                Comments:
        !                  
        !                1) local ghosts comes in two types: ghost that exists
        !                   because two sub domains touch, and ghosts across a 
        !                   periodic boundary. The first is automatically handled
        !                   by the particles - and no copy or/and send is 
        !                   required. The second ghosts could be handled during 
        !                   the calculation of the interactions, but this would 
        !                   require a check for periodicity and no explicit 
        !                   ghosts from the processor itself. However, it seems
        !                   more natural to have ghosts - irrespectively of their
        !                   origin, so ghosts from periodicity (on the same 
        !                   processor) are copied here.
        !
        !                2) is standard procedure
        ! 
        !                3) is handled by copying particles that are adjacent to
        !                   faces on a physical boundary to their image position.
        !
        !                Now, item 2) and 3) can be treated within the same logic
        !                whereas 1) require a bit of thought: to keep the source
        !                compact, we could loop over all processors including
        !                the local one and check for ghosts - the problem with 
        !                this procedure is, that we check for ghosts by comparing
        !                the location of particles within the extended subs 
        !                boundaries (and not their true ghost layer - which 
        !                would result in more IFs than necessary). However, 
        !                because of this, looping over the processor itself 
        !                would find ghosts that are not really ghost - the only 
        !                ghosts that should be found in this step are those due 
        !                to periodicity.  The solution is to store the ghosts as
        !                ighost() with nghost denoting the total number of 
        !                ghosts (excluding those due to periodicity) and 
        !                nghostplus the total number of ghosts. During the loop 
        !                over the processor itself we therefore only need to 
        !                loop from nghost+1,nghostplus.
        !
        ! References  :
        !
        !  Revisions    :
        !------------------------------------------------------------------------
        !  30.11.2009 Xin Bian
        !  When use symmetry inter-communication, make sure that
        !  we create potential ghosts from upper/right physical boundaires of
        !  symmetry, wall_symmetry, Lees-Edwards.
        !
        !  27.07.2009, Xin Bian
        !  make sure the particle lie on a position which is 
        !  exactly ghost_size away from subdomain boundary, 
        !  doesn't need to be a ghost.
        !  When create potential ghost, '=' is excluded from
        !  '<=' or '>='. 
        !  However, we recognize ghost '=' is included into
        !  '<=' or '>=' for cautiousness.  
        !
        ! 
        !  $Log: ppm_map_part_ghost_get.f,v $
        !  Revision 1.16  2006/09/04 18:34:50  pchatela
        !  Fixes and cleanups to make it compile for
        !  release-1-0
        !
        !  Revision 1.15  2006/02/03 09:41:26  ivos
        !  Added the PRELIMINARY ghost_put functionality. Still needs clean-up,
        !  but should work.
        !
        !  Revision 1.14  2005/08/24 15:32:45  hiebers
        !  BUGFIX: nlist2 is also initialized outside loop for the case that
        !  ppm_nsublist = 0
        !
        !  Revision 1.13  2004/11/11 15:26:17  ivos
        !  Moved allocatable work arrays to module.
        !
        !  Revision 1.12  2004/10/01 16:09:06  ivos
        !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
        !
        !  Revision 1.11  2004/08/27 11:40:20  walther
        !  Bug fix in ghosts at right,upper,top boundaries.
        !
        !  Revision 1.10  2004/08/24 12:09:00  walther
        !  Bug fix in allocation of ighost() for periodic systems.
        !
        !  Revision 1.9  2004/08/03 12:46:38  ivos
        !  bugfix: if KIND and ppm_kind were not the same, things went wrong.
        !  Fixed by adding proper IFs and type conversions.
        !
        !  Revision 1.8  2004/07/26 07:42:44  ivos
        !  Changed to single-interface modules and adapted all USE statements.
        !
        !  Revision 1.7  2004/07/15 14:16:57  walther
        !  Bug fix: a ghost could incorrectly be send more than once to a 
        !  processor; introducing the lghost(:) array fixed this problem.
        !
        !  Revision 1.6  2004/07/07 13:09:35  walther
        !  Bug fix: particles in a sub require GE and LT.
        !
        !  Revision 1.5  2004/07/01 13:15:23  walther
        !  Bug fix: allocation of ppm_sendbuffer and ppm_buffer2part.
        !
        !  Revision 1.4  2004/06/03 16:07:22  ivos
        !  Removed debug PRINT statement.
        !
        !  Revision 1.3  2004/06/03 08:49:39  walther
        !  Bug fix: now reallocating the ppm_buffer2part and ppm_sendbuffers/d.
        !  Also skipping the do-loops in step 1 if if we do not have periodicity.
        !
        !  Revision 1.2  2004/05/28 10:32:56  walther
        !  First functional release - without symmetry.
        !
        !  Revision 1.1  2004/02/19 15:56:35  walther
        !  Initial implementation (not complete!).
        !
        !------------------------------------------------------------------------
        !  Parallel Particle Mesh Library (PPM)
        !  Institute of Computational Science
        !  ETH Zentrum, Hirschengraben 84
        !  CH-8092 Zurich, Switzerland
        !-------------------------------------------------------------------------

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_map_part_ghost_get_s(xp,lda,Npart,isymm,ghostsize,info,shear_length)
#elif  __KIND == __DOUBLE_PRECISION
        SUBROUTINE ppm_map_part_ghost_get_d(xp,lda,Npart,isymm,ghostsize,info,shear_length)
#endif        
        !----------------------------------------------------
        !  Modules 
        !----------------------------------------------------
        
        USE ppm_module_data
        USE ppm_module_substart
        USE ppm_module_substop
        USE ppm_module_error
        USE ppm_module_alloc

        IMPLICIT NONE
        
#if    __KIND == __SINGLE_PRECISION  | __KIND_AUX == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#else
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
        
        !----------------------------------------------------
        ! Arguments
        ! shear_length : sheard length in adjacent layer boxes,
        !                used only for Lees-Edwards boundary.
        !----------------------------------------------------
        
        REAL(MK), DIMENSION(:,:), INTENT(IN)    :: xp
        REAL(MK)                , INTENT(IN)    :: ghostsize
        INTEGER                 , INTENT(IN)    :: Npart,lda
        INTEGER                 , INTENT(IN)    :: isymm
        INTEGER                 , INTENT(OUT)   :: info
        REAL(MK),DIMENSION(:,:),OPTIONAL,INTENT(IN) :: shear_length
        
        !----------------------------------------------------
        ! Local variables 
        ! xt                  : position of potential ghosts
        ! xmini, ymini, zmini : inner domain for creating 
        !                       potential ghosts.
        !                       outter domain for recoginzing
        !                       ghosts created by other
        !                       processes.
        ! xmaxi, ymaxi, zmaxi : inner domain for creating 
        !                       potential ghosts.
        !                       outter domain for recoginzing
        !                       ghosts created by other
        !                       processes.
        !----------------------------------------------------
        
        INTEGER, DIMENSION(3) :: ldu
        INTEGER               :: i,j,k
        INTEGER               :: topoid,isub
        INTEGER               :: nlist1,nlist2
        INTEGER               :: nghost,nghostplus
        INTEGER               :: ipart,sendrank,recvrank
        INTEGER               :: iopt,iset,ibuffer
        REAL(MK),DIMENSION(:,:),POINTER :: xt
        INTEGER, DIMENSION(:),POINTER :: tlist
        REAL(MK)              :: xminf,yminf,zminf 
        REAL(MK)              :: xmaxf,ymaxf,zmaxf
        REAL(MK)              :: xmini,ymini,zmini
        REAL(MK)              :: xmaxi,ymaxi,zmaxi
        REAL(MK), DIMENSION(3) :: len_phys
        REAL(MK)              :: t0
        
        !----------------------------------------------------
        ! num_pe    : number of periodic boundaries.
        ! num_sym   : number of symmetry boundaries.
        ! num_wall_sym : number of wall boundaries using
        !                symmetry/mirror particles.
        ! num_le    : number of Lees-Edwards boundaries.
        ! num_extra : num_sym + num_wall_sym + num_le
        !             For symmetry inter-communication
        !             case, we need extra layers
        !             of potential ghosts which are
        !             at physical boundaries.
        ! num_own   : number of boundaries that
        !             may generate ghost particles 
        !             sent to its own process.
        ! lextra    : indicate if we need at each
        !             side of different direction a extra
        !             layer.
        ! min_phys  : minimum boundaries of physics
        !             domain.
        ! max_phys  : maximum boundaries of physics
        !             domain.
        ! min_sub  : minimum boundaries of sub-
        !             domain.
        ! max_sub  : maximum boundaries of sub-
        !             domain. 
        !----------------------------------------------------
        
        INTEGER               :: num_pe
        INTEGER               :: num_sym
        INTEGER               :: num_wall_sym
        INTEGER               :: num_le
        INTEGER               :: num_extra
        LOGICAL, DIMENSION(6) :: lextra
        INTEGER               :: num_own
        REAL(MK), DIMENSION(3):: min_phys
        REAL(MK), DIMENSION(3):: max_phys
        REAL(MK), DIMENSION(3):: min_sub
        REAL(MK), DIMENSION(3):: max_sub
        
        REAL(MK)              :: ppm_machine_zero
        
        
        !----------------------------------------------------
        ! Initialization.
        !----------------------------------------------------
        
        CALL substart('ppm_map_part_ghost_get',t0,info)
        topoid = ppm_topoid
        
        NULLIFY(xt)                
        NULLIFY(tlist)
        
        !----------------------------------------------------
        ! Define machine zero for furture comparison.
        !----------------------------------------------------
        
        ppm_machine_zero = 1.0e-10_MK
        
        IF (MK == ppm_kind_single) THEN
           ppm_machine_zero = ppm_myepss
        ELSE IF (MK == ppm_kind_double ) THEN
           ppm_machine_zero = ppm_myepsd         
        END IF
        
        !----------------------------------------------------
        ! Count number of different boundary conditions.
        !----------------------------------------------------
        
        num_pe       = 0
        num_sym      = 0
        num_wall_sym = 0
        num_le       = 0
        num_extra    = 0
        lextra(:)    = .FALSE.
        num_own      = 0
        
        DO k = 1, 2*ppm_dim
           
           SELECT CASE ( ppm_bcdef(k,topoid) )
              
           CASE ( ppm_param_bcdef_periodic ) 
              
              num_pe    = num_pe + 1
              
           CASE ( ppm_param_bcdef_symmetry )
              
              num_sym   = num_sym + 1
              lextra(k) = .TRUE.
              
           CASE ( ppm_param_bcdef_wall_sym )
              
              num_wall_sym  = num_wall_sym + 1
              lextra(k) = .TRUE.
              
           CASE ( ppm_param_bcdef_LE )
              
              num_le    = num_le + 1
              lextra(k) = .TRUE.
              
           END SELECT
         
        END DO ! k, 1, 2*ppm_dim
        
        num_extra = num_sym + num_wall_sym + num_le
        num_own   = num_pe + num_extra
        
        
        !----------------------------------------------------
        ! Save min and max physics boundaries.
        ! Compute size of computational box on this topology.
        ! Save min and max subdomain boundaries.
        !----------------------------------------------------
        
#if    __KIND == __SINGLE_PRECISION
        
        min_phys(1:ppm_dim) = ppm_min_physs(1:ppm_dim,topoid)
        max_phys(1:ppm_dim) = ppm_max_physs(1:ppm_dim,topoid)
        len_phys(1:ppm_dim) = max_phys(1:ppm_dim) - min_phys(1:ppm_dim)
        
#else
        
        min_phys(1:ppm_dim) = ppm_min_physd(1:ppm_dim,topoid)
        max_phys(1:ppm_dim) = ppm_max_physd(1:ppm_dim,topoid)
        len_phys(1:ppm_dim) = max_phys(1:ppm_dim) - min_phys(1:ppm_dim)
        
#endif         
        
        !----------------------------------------------------
        ! Save the map type for the subsequent calls.
        !----------------------------------------------------
        
        ppm_map_type = ppm_param_map_ghost_get
        
        !----------------------------------------------------
        ! Allocate memory for the list of particles on local
        ! processor that may be ghosts on other processors.
        ! It is at most Npart.
        !
        ! ilist1 having particles' indice for searching.
        ! ilist2 having particles' indice not checked.
        ! ighost saving particle's index when it is potential
        !        ghost.
        !----------------------------------------------------
        
        iopt   = ppm_param_alloc_fit
        ldu(1) = Npart
        
        CALL ppm_alloc(ilist1,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',&
                &        'list1',__LINE__,info)
           GOTO 9999
        ENDIF
        
        CALL ppm_alloc(ilist2,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',&
                &        'list2',__LINE__,info)
           GOTO 9999
        ENDIF
        
        CALL ppm_alloc(ighost,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',&
                &        'ighost',__LINE__,info)
           GOTO 9999
        ENDIF
        
        !----------------------------------------------------
        ! Allocate memory for positions of potential ghosts.
        ! It is ppm_dim dimesnion and at most Npart.
        !----------------------------------------------------
        
        ldu(1) = ppm_dim
        ldu(2) = Npart
        
        CALL ppm_alloc(xt,ldu,iopt,info)
        IF (info.NE.0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',&
                &        'xt',__LINE__,info)
           GOTO 9999
        ENDIF
        
        !----------------------------------------------------
        ! List ilist1() holds the particles we are currently 
        ! considering.
        ! List ilist2() holds the particles that have not yet
        ! been associated with a sub.
        ! ighost() holds the potential  ghosts.
        !----------------------------------------------------
        
        DO i = 1, Npart
           
           ilist1(i) = i
           
        END DO
        
        nlist1  = Npart
        nlist2  = 0
        nghost  = 0
        
        !----------------------------------------------------
        ! Fill the list with particles that are within the
        ! reach of the ghost regions of other subs.
        ! 
        ! We do that by looping over the subs belonging
        ! to this processor and checking if the particles are
        ! within the sub. 
        ! If this is not the case, the particle will never 
        ! be a ghost.
        !
        ! When create these potential ghosts into the list,
        ! we use '<' or '>' without equality, 
        ! since '=ghostsize' would have zero value from kernel,
        ! and therefore no influence to others.
        !
        ! However, later when recognize ghosts created by 
        ! other subs,  we use '<=' or '>=' to be cautiously 
        ! recongizing every ghost particle.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Loop over all subdomains.
        !----------------------------------------------------
        
        DO k = 1, ppm_nsublist(topoid)
           
           !-------------------------------------------------
           ! Initialize the second list counter to zero.
           !-------------------------------------------------
           
           nlist2 = 0
           
           !-------------------------------------------------
           ! Get the (global) id of this sub.
           !-------------------------------------------------
           
           isub = ppm_isublist(k,topoid)
           
           !-------------------------------------------------
           ! Store the extend of the sub in minf and maxf.
           !-------------------------------------------------

#if    __KIND == __SINGLE_PRECISION
                      
           xminf = ppm_min_subs(1,isub,topoid)
           xmaxf = ppm_max_subs(1,isub,topoid)
           
           yminf = ppm_min_subs(2,isub,topoid)
           ymaxf = ppm_max_subs(2,isub,topoid)
           
           IF ( ppm_dim == 3 ) THEN
              zminf = ppm_min_subs(3,isub,topoid)
              zmaxf = ppm_max_subs(3,isub,topoid)
           END IF
        
#else
           
           xminf = ppm_min_subd(1,isub,topoid)
           xmaxf = ppm_max_subd(1,isub,topoid)
           
           yminf = ppm_min_subd(2,isub,topoid)
           ymaxf = ppm_max_subd(2,isub,topoid)
           
           IF ( ppm_dim == 3 ) THEN
              zminf = ppm_min_subd(3,isub,topoid)
              zmaxf = ppm_max_subd(3,isub,topoid)
           END IF

#endif   
           
           !-------------------------------------------------
           ! Compute the size of the inner region.
           !-------------------------------------------------
         
           IF ( isymm > 0 ) THEN
              
              !----------------------------------------------
              ! if we use symmetry, for periodic and solid 
              ! wall boundaries, the upper/right part of our 
              ! sub will have ghosts and therefore particles
              ! at the upper/right part of the sub cannot be
              ! ghosts on other processors. Thus the ghosts
              ! must be found at the lower/left of the sub.
              !
              ! But, for symmetry, wall_symmetry or 
              ! Lees-Edwards boundary, the particles along
              ! lower/left part of physics boundary can
              ! have ghosts also, therefore,
              ! upper/right part of physical boundary
              ! can be potential ghosts also. 
              ! Therefore, we create potential ghosts along 
              ! these upper/right boundaries.
              !----------------------------------------------
              
              xmini = xminf + ghostsize
              xmaxi = xmaxf
              
              ymini = yminf + ghostsize
              ymaxi = ymaxf
              
              IF ( ppm_dim == 3 ) THEN
                 
                 zmini = zminf + ghostsize 
                 zmaxi = zmaxf
                 
              END IF
              
              !----------------------------------------------
              ! If this sub-domain is at physics boundary and
              ! we need extra layer of ghosts in this boundary,
              ! then we extend inner region of finding
              ! potential ghosts.
              !----------------------------------------------
              
              IF ( ABS( xmaxi - max_phys(1) ) < ppm_machine_zero .AND. &
                   lextra(2) ) THEN
                 
                 xmaxi = xmaxi - ghostsize
                 
              END IF
              
              IF ( ABS( ymaxi - max_phys(2) ) < ppm_machine_zero .AND. &
                   lextra(4) ) THEN
                 
                 ymaxi = ymaxi - ghostsize
                 
              END IF
              
              IF ( ppm_dim == 3 .AND. &
                   ABS( zmaxi - max_phys(3) ) < ppm_machine_zero .AND. &
                   lextra(6) ) THEN
                 
                 zmaxi = zmaxi - ghostsize
                 
              END IF
              
           ELSE
              
              !----------------------------------------------
              ! If we do not use symmetry,
              ! the particles along the boundary of the sub
              ! will be potential ghosts on other processes.
              !----------------------------------------------
              
              xmini = xminf + ghostsize
              xmaxi = xmaxf - ghostsize
              
              ymini = yminf + ghostsize
              ymaxi = ymaxf - ghostsize
              
              IF ( ppm_dim ==3  ) THEN
                 
                 zmini = zminf + ghostsize
                 zmaxi = zmaxf - ghostsize
                 
              END IF
              
           END IF
           
           !-------------------------------------------------
           ! loop over the remaining particles, and
           !  record those lying inside ghost regions.
           !-------------------------------------------------
           
           IF ( ppm_dim == 2 ) THEN
              
              DO j = 1, nlist1
                 
                 !-------------------------------------------
                 ! get the particle index.
                 !-------------------------------------------
                 
                 ipart = ilist1(j)
                 
                 !-------------------------------------------
                 ! check if the particle belongs to this sub,
                 ! including '=' for inequality check.
                 !-------------------------------------------
                 
                 IF ( xp(1,ipart) >= xminf .AND. &
                      xp(1,ipart) <= xmaxf .AND. &
                      xp(2,ipart) >= yminf .AND. &
                      xp(2,ipart) <= ymaxf ) THEN
                    
                    !----------------------------------------
                    ! Check if particles are in potential 
                    ! ghost layer, without '=' for inequality
                    ! check.
                    !----------------------------------------
                    
                    IF ( xp(1,ipart) < xmini .OR. &
                         xp(1,ipart) > xmaxi .OR. &
                         xp(2,ipart) < ymini .OR. &
                         xp(2,ipart) > ymaxi ) THEN
                       
                       !-------------------------------------
                       ! The particle may be a ghost
                       ! somewhere, record it.
                       !-------------------------------------
                       
                       nghost         = nghost + 1
                       ighost(nghost) = ipart
                       xt(1,nghost)   = xp(1,ipart)
                       xt(2,nghost)   = xp(2,ipart)
                       
                    END IF
                    
                 ELSE    
                    
                    !----------------------------------------
                    ! If not on this sub,
                    ! we need to consider it further.
                    !----------------------------------------
                    
                    nlist2         = nlist2 + 1
                    ilist2(nlist2) = ipart
                    
                 END IF
                 
              END DO
              
           ELSE ! ppm_dim == 3D
              
              DO j = 1 , nlist1
                 
                 !-------------------------------------------
                 ! get the particle index.
                 !-------------------------------------------
                 
                 ipart = ilist1(j)
                 
                 !-------------------------------------------
                 ! check if the particle belongs to this sub,
                 ! including '=' for inequality check.
                 !-------------------------------------------
                 
                 IF ( xp(1,ipart) >= xminf .AND. &
                      xp(1,ipart) <= xmaxf .AND. &
                      xp(2,ipart) >= yminf .AND. &
                      xp(2,ipart) <= ymaxf .AND. &
                      xp(3,ipart) >= zminf .AND. &
                      xp(3,ipart) <= zmaxf ) THEN
                    
                    !----------------------------------------
                    ! check if the particle is in potential 
                    ! ghost layer, without '=' for inequality
                    ! check.
                    !----------------------------------------
                    
                    IF ( xp(1,ipart) < xmini .OR. &
                         xp(1,ipart) > xmaxi .OR. &
                         xp(2,ipart) < ymini .OR. &
                         xp(2,ipart) > ymaxi .OR. &
                         xp(3,ipart) < zmini .OR. &
                         xp(3,ipart) > zmaxi ) THEN
                       
                       !-------------------------------------
                       ! The particle may be a ghost
                       ! somewhere, record its index and
                       ! position.
                       !-------------------------------------
                       
                       nghost         = nghost + 1   
                       ighost(nghost) = ipart
                       xt(1,nghost)   = xp(1,ipart)
                       xt(2,nghost)   = xp(2,ipart)
                       xt(3,nghost)   = xp(3,ipart)
                       
                    END IF
                    
                 ELSE    
                    
                    !----------------------------------------
                    ! If not on this sub we need to save
                    ! its index and consider further.
                    !----------------------------------------
                    
                    nlist2         = nlist2 + 1
                    ilist2(nlist2) = ipart
                    
                 END IF
                 
              END DO ! j = 1, ilist1
              
           END IF ! ppm_dim == 3
           
           !-------------------------------------------------
           ! swap the lists.
           ! Changed to use pointer operations, which
           ! is supposed to be more efficient than
           ! array copying.
           !-------------------------------------------------
           
           tlist  => ilist1
           ilist1 => ilist2
           ilist2 => tlist
           
           nlist1 = nlist2
           
        END DO ! k = 1,  ppm_nsublist()
        
        
        !----------------------------------------------------
        ! At the end, the nlist2 should be zero,
        ! since particles on this process should belong one
        ! of the subdomains on this process.
        !----------------------------------------------------
        
        IF ( nlist2 .NE.0 ) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_part_unass,'ppm_map_part_ghost_get',&
                'nlist2 > 0',__LINE__,info)
           GOTO 9999
        END IF
        
        
        !----------------------------------------------------
        ! Ok, we now have a list of potential ghosts. 
        ! From these we extract/add their periodic images,
        ! or symmetry boundary particles, or wall_symmetry 
        ! boundary particles, 
        ! or Lees-Edwards boundary particles.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Initialize the total number potential ghosts,
        ! which doesn't include extra ghosts yet.
        !----------------------------------------------------
        
        nghostplus = nghost
        

        !----------------------------------------------------
        ! Lees-Edwards boundary.
        !----------------------------------------------------
        
        IF ( ppm_bcdef(1,topoid) == ppm_param_bcdef_LE .AND. &
             PRESENT(shear_length) ) THEN
           
#include "ppm_bcdef_le_x.inc"
           
        END IF
        
        IF ( ppm_bcdef(3,topoid) == ppm_param_bcdef_LE .AND. &
             PRESENT(shear_length) ) THEN
           
#include "ppm_bcdef_le_y.inc"
           
        END IF
        
        IF ( ppm_dim == 3 ) THEN
           
           IF ( ppm_bcdef(5,topoid) == ppm_param_bcdef_LE .AND. &
                PRESENT(shear_length) ) THEN
              
#include "ppm_bcdef_le_z.inc"
              
           END IF
           
        END IF
        
        
        !----------------------------------------------------
        ! Wall boundary outside of computational domain 
        ! using symmetry/mirror particles.
        !----------------------------------------------------
        
        IF ( ppm_bcdef(1,topoid) == ppm_param_bcdef_wall_sym ) THEN
           
#include "ppm_bcdef_symmetry_lower_x.inc"
           
        END IF
        
        IF ( ppm_bcdef(2,topoid) == ppm_param_bcdef_wall_sym ) THEN
           
#include "ppm_bcdef_symmetry_upper_x.inc"
           
        END IF
        
        IF ( ppm_bcdef(3,topoid) == ppm_param_bcdef_wall_sym ) THEN
           
#include "ppm_bcdef_symmetry_lower_y.inc"
           
        END IF
        
        IF ( ppm_bcdef(4,topoid) == ppm_param_bcdef_wall_sym ) THEN
           
#include "ppm_bcdef_symmetry_upper_y.inc"
           
        ENDIF
        
        IF ( ppm_dim == 3 ) THEN
           
           IF ( ppm_bcdef(5,topoid) == ppm_param_bcdef_wall_sym ) THEN
              
#include "ppm_bcdef_symmetry_lower_z.inc"
              
           END IF
           
           IF ( ppm_bcdef(6,topoid) == ppm_param_bcdef_wall_sym ) THEN
              
#include "ppm_bcdef_symmetry_upper_z.inc"
              
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Symmetry boundary.
        !----------------------------------------------------
        
        IF ( ppm_bcdef(1,topoid) == ppm_param_bcdef_symmetry ) THEN
           
#include "ppm_bcdef_symmetry_lower_x.inc"
           
        END IF
        
        IF ( ppm_bcdef(2,topoid) == ppm_param_bcdef_symmetry ) THEN
           
#include "ppm_bcdef_symmetry_upper_x.inc"
           
        END IF
        
        IF ( ppm_bcdef(3,topoid) == ppm_param_bcdef_symmetry ) THEN
           
#include "ppm_bcdef_symmetry_lower_y.inc"
           
        END IF
        
        IF ( ppm_bcdef(4,topoid) == ppm_param_bcdef_symmetry ) THEN
           
#include "ppm_bcdef_symmetry_upper_y.inc"
           
        END IF
        
        IF ( ppm_dim == 3 ) THEN
           
           IF ( ppm_bcdef(5,topoid) == ppm_param_bcdef_symmetry ) THEN
              
#include "ppm_bcdef_symmetry_lower_z.inc"
              
           END IF
           
           IF ( ppm_bcdef(6,topoid) == ppm_param_bcdef_symmetry ) THEN
              
#include "ppm_bcdef_symmetry_upper_z.inc"
              
           END IF
        END IF
        
        !----------------------------------------------------
        ! Periodic boundary.
        !----------------------------------------------------
        
        IF ( ppm_bcdef(1,topoid) == ppm_param_bcdef_periodic ) THEN
           
#include "ppm_bcdef_periodic_x.inc"
           
        END IF
        
        IF (ppm_bcdef(3,topoid) == ppm_param_bcdef_periodic ) THEN
           
#include "ppm_bcdef_periodic_y.inc"
           
        END IF
        
        IF ( ppm_dim == 3 ) THEN
           
           IF ( ppm_bcdef(5,topoid) == ppm_param_bcdef_periodic ) THEN
              
#include "ppm_bcdef_periodic_z.inc"
              
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Now we have a list of all potential ghosts to search,
        ! allocate memory for the lghosts, which is a logical
        ! array indicating if particle being taken.
        !----------------------------------------------------
        
        iopt   = ppm_param_alloc_fit
        ldu(1) = nghostplus
        
        CALL ppm_alloc(lghost,ldu,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
                'logical list: lghost',__LINE__,info)
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Allocate memory for the global send/recv lists,
        ! ppm_ncommseq has the total number of ranks.
        !----------------------------------------------------
        
        ldu(1) = ppm_ncommseq(topoid)
        
        CALL ppm_alloc(ppm_isendlist,ldu,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
                'global send rank list PPM_ISENDLIST',__LINE__,info)
           GOTO 9999
        END IF
        
        CALL ppm_alloc(ppm_irecvlist,ldu,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
                'global recv rank list PPM_IRECVLIST',__LINE__,info)
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Allocate memory for pointers to the particles that 
        ! will be sent to and received from a process.
        !
        ! The ppm_recvbuffer is NOT used in this routine, 
        ! but initialized from the precv() array in the routine
        ! ppm_map_part_send().
        !
        ! ppm_ncommseq() has the number of total ranks.
        !----------------------------------------------------
        
        ldu(1) = ppm_ncommseq(topoid) + 1
        
        CALL ppm_alloc(ppm_psendbuffer,ldu,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
                'global send buffer pointer PPM_PSENDBUFFER', &
                __LINE__,info)
           GOTO 9999
        ENDIF
        
        !----------------------------------------------------
        ! Step 1:
        ! First we find ghosts, due to periodic, symmetry, 
        ! wall_symmetry, Lees-Edwards boundary, that can be on
        ! the local process.
        !----------------------------------------------------
        
        ppm_psendbuffer(1) = 1
        ppm_nsendlist      = 0
        ppm_nrecvlist      = 0
        iset               = 0
        ibuffer            = 0
        k                  = 1
        
        !----------------------------------------------------
        ! Get the rank of the process.
        ! ppm_icommseq(1,topoid) is the rank of current
        ! process.
        !----------------------------------------------------
        
        sendrank = ppm_icommseq(k,topoid) ! should be ppm_rank.
        recvrank = sendrank
        
        !----------------------------------------------------
        ! Store the process to which we will send.
        !----------------------------------------------------
        
        ppm_nsendlist                = ppm_nsendlist + 1
        ppm_isendlist(ppm_nsendlist) = sendrank
        
        !----------------------------------------------------
        ! Store the processor from which we will recv.
        !----------------------------------------------------
        
        ppm_nrecvlist                = ppm_nrecvlist + 1
        ppm_irecvlist(ppm_nrecvlist) = recvrank
        
        !----------------------------------------------------
        ! Store the number of buffer entries.
        ! (this is the first)
        !----------------------------------------------------
        
        ppm_buffer_set = 1
        
        !----------------------------------------------------
        ! Allocate memory for the field registers that holds 
        ! the dimension and type of the data
        !----------------------------------------------------
        
        iopt   = ppm_param_alloc_fit
        ldu(1) = ppm_buffer_set
        CALL ppm_alloc(ppm_buffer_dim ,ldu,iopt,info)
        
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
                'buffer dimensions PPM_BUFFER_DIM',&
                __LINE__,info)
           GOTO 9999
        END IF
        
        CALL ppm_alloc(ppm_buffer_type,ldu,iopt,info)
        
        IF (info .NE. 0) THEN
           info = ppm_error_fatal
           CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',&
                'buffer types PPM_BUFFER_TYPE',&
                __LINE__,info)
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Store the dimension and type.
        !----------------------------------------------------
        
        ppm_buffer_dim(ppm_buffer_set)  = ppm_dim
#if    __KIND == __SINGLE_PRECISION
        ppm_buffer_type(ppm_buffer_set) = ppm_kind_single
#else
        ppm_buffer_type(ppm_buffer_set) = ppm_kind_double
#endif
        
        !----------------------------------------------------
        ! Well we only do the stuff for real if we have 
        ! periodic, symmetry, wall_symmetry, or
        ! Lees-Edwards boundary.
        !
        ! (Re)allocate memory for send buffer. The required 
        ! size of these arrays is not easy to compute.
        !
        ! The first step (step 1) require at most 
        ! ppm_nsublist(topoid)*(nghostplus - nghost)*ppm_dim.
        !
        ! The arrays are resized further during step 2 and 3.
        !----------------------------------------------------
        
        IF ( num_own > 0 ) THEN
           
           !-------------------------------------------------
           ! We can grow or fit the arrays, a matter of taste.
           !-------------------------------------------------
           
           iopt   = ppm_param_alloc_grow
           ldu(1) = ppm_dim * (nghostplus - nghost) * &
                ppm_nsublist(topoid)
           
           !-------------------------------------------------
           ! First allocate the sendbuffer.
           !-------------------------------------------------
           
           IF ( ppm_kind == ppm_kind_double ) THEN
              
              CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
              
           ELSE
              
              CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
              
           END IF
           
           IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
                   'global send buffer PPM_SENDBUFFER', &
                   __LINE__,info)
              GOTO 9999
           ENDIF
           
           !-------------------------------------------------
           ! then allocate the index list: buffer2part.
           !-------------------------------------------------
           
           ldu(1) = (nghostplus - nghost) * ppm_nsublist(topoid)
           CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
           IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
                   'buffer-to-particles map PPM_BUFFER2PART', &
                   __LINE__,info)
              GOTO 9999
           ENDIF
           
           !-------------------------------------------------
           ! flag all ghosts as not yet taken.
           !-------------------------------------------------
           
           lghost(:) = .TRUE.
           
           !-------------------------------------------------
           ! loop over the subs on the local process.
           !-------------------------------------------------
           
           DO j = 1, ppm_nsublist(topoid)
              
              !----------------------------------------------
              ! Get the global ID of the sub.
              !----------------------------------------------
              
              isub = ppm_isublist(j,topoid) 
              
              !----------------------------------------------
              ! Define the extended resize of this sub 
              !----------------------------------------------

#if    __KIND == __SINGLE_PRECISION
              
              xminf = ppm_min_subs(1,isub,topoid)
              xmaxf = ppm_max_subs(1,isub,topoid)
              
              yminf = ppm_min_subs(2,isub,topoid)
              ymaxf = ppm_max_subs(2,isub,topoid)
              
              IF ( ppm_dim == 3 ) THEN
                 zminf = ppm_min_subs(3,isub,topoid)
                 zmaxf = ppm_max_subs(3,isub,topoid)
              END IF
              
#else
              
              xminf = ppm_min_subd(1,isub,topoid)
              xmaxf = ppm_max_subd(1,isub,topoid)
              
              yminf = ppm_min_subd(2,isub,topoid)
              ymaxf = ppm_max_subd(2,isub,topoid)
              
              IF ( ppm_dim == 3 ) THEN
                 zminf = ppm_min_subd(3,isub,topoid)
                 zmaxf = ppm_max_subd(3,isub,topoid)
              END IF
              
#endif
              
              IF ( isymm > 0 ) THEN
                 
                 !-------------------------------------------
                 ! if we use symmetry ghosts will only be 
                 ! present at the upper/right part of the sub.
                 !-------------------------------------------
                 
                 xmini = xminf
                 xmaxi = xmaxf + ghostsize
                 
                 ymini = yminf
                 ymaxi = ymaxf + ghostsize
                 
                 IF ( ppm_dim == 3 ) THEN
                    
                    zmini = zminf
                    zmaxi = zmaxf + ghostsize
                    
                 END IF
                 
                 !-------------------------------------------
                 ! If this sub domain is at physics boundary,
                 ! and we need extra layer of ghosts in this
                 ! boundary. 
                 !-------------------------------------------
                 
                 IF ( ABS( xmini - min_phys(1) ) < ppm_machine_zero .AND. &
                      lextra(1) ) THEN
                    xmini = xmini - ghostsize
                 END IF
                 
                 IF ( ABS( ymini - min_phys(2) ) < ppm_machine_zero .AND. &
                      lextra(3) ) THEN
                    ymini = ymini - ghostsize
                 END IF
                 
                 IF  ( ppm_dim == 3 .AND. &
                      ABS( zmini - min_phys(3) ) < ppm_machine_zero .AND. &
                      lextra(5) ) THEN
                    zmini = zmini - ghostsize
                 END IF
                 
              ELSE
                 
                 !-------------------------------------------
                 ! if we do not use symmetry, we have ghost 
                 ! all around sub-domain's boundaries.
                 !-------------------------------------------
                 
                 xmini = xminf - ghostsize
                 xmaxi = xmaxf + ghostsize
                 
                 ymini = yminf - ghostsize
                 ymaxi = ymaxf + ghostsize
                 
                 IF ( ppm_dim == 3 ) THEN
                    zmini = zminf - ghostsize
                    zmaxi = zmaxf + ghostsize
                 END IF
                 
              END IF
              
              !----------------------------------------------
              ! Considering more than one subdomains in each
              ! process, reallocation of memory has to be done.
              !
              ! Reallocate the arrays: 
              ! ppm_sendbuffers/d and ppm_buffer2part. 
              !
              ! The required size is current size minus 
              ! current use plus maximum increment
              ! (nghostplus-nghost).
              !
              ! We perform a grow_preserve to preserve the 
              ! contents but only to reallocate if an 
              ! increased size is required.
              !----------------------------------------------
              
              iopt = ppm_param_alloc_grow_preserve
              
              IF ( ppm_kind == ppm_kind_double ) THEN
                 
                 IF ( ASSOCIATED(ppm_sendbufferd) ) THEN
                    ldu(1) = SIZE(ppm_sendbufferd) + &
                         (nghostplus - nghost - ibuffer)*ppm_dim
                 ELSE
                    ldu(1) = (nghostplus - ibuffer)*ppm_dim
                 END IF
                 
                 CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
                 
              ELSE
                 
                 IF ( ASSOCIATED(ppm_sendbuffers) ) THEN
                    ldu(1) = SIZE(ppm_sendbuffers) + &
                         (nghostplus - nghost - ibuffer)*ppm_dim
                 ELSE
                    ldu(1) = (nghostplus - nghost - ibuffer)*ppm_dim 
                 END IF
                 
                 CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
                 
              END IF
              
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
                      'global send buffer PPM_SENDBUFFER',&
                      __LINE__,info)
                 GOTO 9999
              ENDIF
              
              IF ( ASSOCIATED(ppm_buffer2part) ) THEN
                 ldu(1) = SIZE(ppm_buffer2part) + nghostplus - nghost - iset
              ELSE
                 !---------------------------------
                 ! Changed by Xin 
                 ! ldu(1) = nghostplus - nghost - iset
                 !---------------------------------
                 ldu(1) = nghostplus  - iset
              END IF
              
              CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
              
              IF (info .NE. 0) THEN
                 info = ppm_error_fatal
                 CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',&
                      'buffer-to-particles map PPM_BUFFER2PART',&
                      __LINE__,info)
                 GOTO 9999
              END IF
              
              !----------------------------------------------
              ! loop over potential ghost particles 
              ! due to periodicity, symmetry, wall_symmetry,
              ! Lees-Edwards boundary.
              !----------------------------------------------
              
              IF ( ppm_dim == 2 ) THEN
                 
                 !-------------------------------------------
                 ! Two dimensions.
                 !-------------------------------------------
                 
                 DO i = nghost+1, nghostplus
                    
                    !----------------------------------------
                    ! check if it is inside the ghost region.
                    !----------------------------------------
                    
                    IF ( xt(1,i) >= xmini .AND. &
                         xt(1,i) <= xmaxi .AND. &
                         xt(2,i) >= ymini .AND. &
                         xt(2,i) <= ymaxi .AND. &
                         lghost(i) ) THEN
                       
                       !-------------------------------------
                       ! mark the ghost as taken.
                       !-------------------------------------
                       
                       lghost(i) = .FALSE.
                       
                       !-------------------------------------
                       ! found one, increment buffer counter.
                       !-------------------------------------
                       
                       iset  = iset + 1
                       
                       !-------------------------------------
                       ! store the ID of the particle.
                       !-------------------------------------
                       
                       ppm_buffer2part(iset) = ighost(i)
                       
                       !-------------------------------------
                       ! store the particle's position.
                       !-------------------------------------
                       
                       IF ( ppm_kind == ppm_kind_double ) THEN
                          
                          ibuffer = ibuffer + 1 
                          ppm_sendbufferd(ibuffer) = xt(1,i)
                          ibuffer = ibuffer + 1 
                          ppm_sendbufferd(ibuffer) = xt(2,i)
                          
                       ELSE
                          
                          ibuffer = ibuffer + 1 
                          ppm_sendbuffers(ibuffer) = xt(1,i)
                          ibuffer = ibuffer + 1 
                          ppm_sendbuffers(ibuffer) = xt(2,i)
                          
                       END IF
                       
                    END IF ! xt 
                    
                 END DO ! i = nghost+1, nghostplus
                 
              ELSE ! 3D
                 
                 !-------------------------------------------
                 ! Three dimensions.
                 !-------------------------------------------
                 
                 DO i = nghost+1, nghostplus
                    
                    !----------------------------------------
                    ! check if it is inside the ghost region .
                    !----------------------------------------
                    
                    IF ( xt(1,i) >= xmini .AND. &
                         xt(1,i) <= xmaxi .AND. &
                         xt(2,i) >= ymini .AND. &
                         xt(2,i) <= ymaxi .AND. & 
                         xt(3,i) >= zmini .AND. &
                         xt(3,i) <= zmaxi .AND. &
                         lghost(i)) THEN
                       
                       !-------------------------------------
                       ! mark the ghost as taken.
                       !-------------------------------------
                       
                       lghost(i) = .FALSE.
                       
                       !-------------------------------------
                       ! found one, increment buffer counter.
                       !-------------------------------------
                       
                       iset = iset + 1
                       
                       !-------------------------------------
                       ! store the ID of the particle.
                       !-------------------------------------
                       
                       ppm_buffer2part(iset) = ighost(i)
                       
                       !-------------------------------------
                       ! store the particle's position.
                       !-------------------------------------
                       
                       IF ( ppm_kind == ppm_kind_double ) THEN
                          
                          ibuffer = ibuffer + 1 
                          ppm_sendbufferd(ibuffer) = xt(1,i)
                          ibuffer = ibuffer + 1 
                          ppm_sendbufferd(ibuffer) = xt(2,i)
                          ibuffer = ibuffer + 1 
                          ppm_sendbufferd(ibuffer) = xt(3,i)
                          
                       ELSE
                          
                          ibuffer = ibuffer + 1 
                          ppm_sendbuffers(ibuffer) = xt(1,i)
                          ibuffer = ibuffer + 1 
                          ppm_sendbuffers(ibuffer) = xt(2,i)
                          ibuffer = ibuffer + 1 
                          ppm_sendbuffers(ibuffer) = xt(3,i)                          
                          
                       END IF
                       
                    END IF
                    
                 END DO ! i = nghost +1, nghostplus
                 
              END IF ! 2/3 dimension 
              
           END DO ! j ( end of loop over subs on local processor)
           
        END IF ! num_own .GT. 0 ( if there are periodic, symmetry etc. ghosts)
        
        
        !----------------------------------------------------
        ! Update the buffer pointer (i.e., the current iset
        ! or the starting index of particles that will be sent
        ! to the k-th entry in the icommseq list.
        !----------------------------------------------------
        
        ppm_psendbuffer(k+1) = iset + 1
        
        !----------------------------------------------------
        ! Step 2 and 3:
        !
        ! Notice the ppm_icommseq(1,topoid) points to the 
        ! local processor, and we already considered this 
        ! on step 1, so we now start at k = 2
        !----------------------------------------------------
        
        DO k = 2, ppm_ncommseq(topoid)
           
           !-------------------------------------------------
           ! for each processor in the sendrank list we need 
           ! to send the ghosts only once, so we flag the 
           ! lghost as true initially.
           !-------------------------------------------------
           
           lghost(:) = .TRUE.
           
           !-------------------------------------------------
           ! Get the rank of the processor
           !-------------------------------------------------
           
           sendrank = ppm_icommseq(k,topoid)
           recvrank = sendrank
           
           !-------------------------------------------------
           ! Store the process to which we will send.
           !-------------------------------------------------
           
           ppm_nsendlist                = ppm_nsendlist + 1
           ppm_isendlist(ppm_nsendlist) = sendrank
           
           !-------------------------------------------------
           ! Store the process from which we will recv.
           !-------------------------------------------------
           
           ppm_nrecvlist                = ppm_nrecvlist + 1
           ppm_irecvlist(ppm_nrecvlist) = recvrank
           
           !-------------------------------------------------
           ! only consider non-negative sendranks.
           !-------------------------------------------------
           
           IF ( sendrank >= 0 ) THEN
              
              !----------------------------------------------
              ! loop over all subs in the current topology,
              ! this is a bit tedious and perhaps inefficient,
              ! but we do not have the subs belonging to all 
              ! processors stored for each topoid. sorry ! 
              ! so we have to perform the calculations each time.
              !
              ! ppm_nsubs(topoid) is the number of all subdomains
              ! on all processes.
              !----------------------------------------------
              
              DO j = 1, ppm_nsubs(topoid) 
                 
                 !-------------------------------------------
                 ! consider only those beloging to the rank.
                 !-------------------------------------------
                 
                 IF ( ppm_subs2proc(j,topoid) == sendrank ) THEN
                    
                    !----------------------------------------
                    ! Define the extended resize of this sub.
                    !----------------------------------------
                    
#if    __KIND == __SINGLE_PRECISION
                    
                    xminf = ppm_min_subs(1,j,topoid)
                    xmaxf = ppm_max_subs(1,j,topoid)
                    
                    yminf = ppm_min_subs(2,j,topoid)
                    ymaxf = ppm_max_subs(2,j,topoid)
                    
                    IF ( ppm_dim == 3 ) THEN
                       zminf = ppm_min_subs(3,j,topoid)
                       zmaxf = ppm_max_subs(3,j,topoid)
                    END IF
                    
#else
                    
                    xminf = ppm_min_subd(1,j,topoid)
                    xmaxf = ppm_max_subd(1,j,topoid)
                    
                    yminf = ppm_min_subd(2,j,topoid)
                    ymaxf = ppm_max_subd(2,j,topoid)
                    
                    IF ( ppm_dim == 3 ) THEN
                       zminf = ppm_min_subd(3,j,topoid)
                       zmaxf = ppm_max_subd(3,j,topoid)
                    END IF
                    
#endif
                    
                    IF ( isymm > 0 ) THEN
                       
                       !-------------------------------------
                       ! if we use symmetry, ghosts will only
                       ! be present at upper/right part of sub.
                       !-------------------------------------
                       
                       xmini = xminf
                       xmaxi = xmaxf + ghostsize
                       
                       ymini = yminf
                       ymaxi = ymaxf + ghostsize
                       
                       IF ( ppm_dim == 3 ) THEN
                          
                          zmini = zminf
                          zmaxi = zmaxf + ghostsize
                          
                       END IF
                       
                       !-------------------------------------
                       ! If this sub domain is at physics 
                       ! boundary and we need extra layer of 
                       ! ghosts in this boundary. 
                       !-------------------------------------
                       
                       IF ( ABS( xmini - min_phys(1) ) < ppm_machine_zero .AND. &
                            lextra(1) ) THEN
                          
                          xmini = xmini - ghostsize
                          
                       END IF
                       
                       IF ( ABS( ymini - min_phys(2) ) < ppm_machine_zero .AND. &
                            lextra(3) ) THEN
                          
                          ymini = ymini - ghostsize
                          
                       END IF
                       
                       IF  (ppm_dim == 3 .AND. &
                            ABS( zmini - min_phys(3) ) < ppm_machine_zero .AND. &
                            lextra(5) ) THEN
                          
                          zmini = zmini - ghostsize
                          
                       END IF
                       
                    ELSE
                       
                       !-------------------------------------
                       ! if we do not use symmetry, 
                       ! we have ghost all around.
                       !-------------------------------------
                       
                       xmini = xminf - ghostsize
                       xmaxi = xmaxf + ghostsize
                       
                       ymini = yminf - ghostsize
                       ymaxi = ymaxf + ghostsize
                       
                       IF ( ppm_dim == 3 ) THEN

                          zmini = zminf - ghostsize
                          zmaxi = zmaxf + ghostsize
                          
                       END IF
                       
                    END IF
                    
                    !----------------------------------------
                    ! Reallocate to make sure we have enough 
                    ! memory in ppm_buffer2part and 
                    ! ppm_sendbuffers/d.
                    !----------------------------------------
                    
                    iopt   = ppm_param_alloc_grow_preserve
                    ldu(1) = ibuffer + nghostplus*ppm_dim
                    
                    IF ( ppm_kind == ppm_kind_double ) THEN
                       CALL ppm_alloc(ppm_sendbufferd,ldu,iopt,info)
                    ELSE
                       CALL ppm_alloc(ppm_sendbuffers,ldu,iopt,info)
                    END IF
                    IF (info .NE. 0) THEN
                       info = ppm_error_fatal
                       CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',&
                            'global send buffer PPM_SENDBUFFER',__LINE__,info)
                       GOTO 9999
                    ENDIF
                    
                    ldu(1) = iset + nghostplus
                    CALL ppm_alloc(ppm_buffer2part,ldu,iopt,info)
                    IF (info .NE. 0) THEN
                       info = ppm_error_fatal
                       CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get',&
                            'buffer2particles map PPM_BUFFER2PART',__LINE__,info)
                       GOTO 9999
                    END IF
                    
                    !----------------------------------------
                    ! loop over the potential ghost particles. 
                    !----------------------------------------
                    
                    IF ( ppm_dim == 2 ) THEN
                       
                       !-------------------------------------
                       ! Two dimensions.
                       !-------------------------------------
                       
                       DO i = 1, nghostplus
                          
                          !----------------------------------
                          ! check if it is inside ghost region. 
                          !----------------------------------
                          
                          IF ( xt(1,i) >= xmini .AND. &
                               xt(1,i) <= xmaxi .AND. &
                               xt(2,i) >= ymini .AND. &
                               xt(2,i) <= ymaxi .AND. &
                               lghost(i)) THEN
                             
                             !-------------------------------
                             ! Mark the ghost as taken for 
                             ! this sendrank.
                             !-------------------------------
                             
                             lghost(i) = .FALSE.
                             
                             !-------------------------------
                             ! found one, increment buffer counter.
                             !-------------------------------
                             
                             iset  = iset + 1

                             !-------------------------------
                             ! store the ID of the particle.
                             !-------------------------------
                             
                             ppm_buffer2part(iset) = ighost(i)
                             
                             !-------------------------------
                             ! store the particle's position.
                             !---------------------------------
                             
                             IF ( ppm_kind == ppm_kind_double ) THEN
                                
                                ibuffer = ibuffer + 1 
                                ppm_sendbufferd(ibuffer) = xt(1,i)
                                ibuffer = ibuffer + 1 
                                ppm_sendbufferd(ibuffer) = xt(2,i)
                                
                             ELSE
                                
                                ibuffer = ibuffer + 1 
                                ppm_sendbuffers(ibuffer) = xt(1,i)
                                ibuffer = ibuffer + 1 
                                ppm_sendbuffers(ibuffer) = xt(2,i)
                                
                             END IF ! ppm_kind
                             
                          END IF ! xt
                          
                       END DO ! i = 1, nghostplus
                       
                    ELSE
                       
                       !-------------------------------------
                       ! Three dimensions.
                       !-------------------------------------
                       
                       DO i = 1, nghostplus
                          
                          !----------------------------------
                          ! check if it is inside ghost region.
                          !----------------------------------
                          
                          IF ( xt(1,i) >= xmini .AND. &
                               xt(1,i) <= xmaxi .AND. &
                               xt(2,i) >= ymini .AND. &
                               xt(2,i) <= ymaxi .AND. & 
                               xt(3,i) >= zmini .AND. &
                               xt(3,i) <= zmaxi .AND. &
                               lghost(i) ) THEN
                             
                             !-------------------------------
                             ! Mark ghost as taken for this
                             ! sendrank.
                             !-------------------------------
                             
                             lghost(i) = .FALSE.
                             
                             !--------------------------------
                             ! found one, increment buffer 
                             ! counter.
                             !-------------------------------
                             
                             iset   = iset + 1
                             
                             !-------------------------------
                             ! store the ID of the particle.
                             !-------------------------------
                             
                             ppm_buffer2part(iset) = ighost(i)
                             
                             !-------------------------------
                             ! store the particle's position.
                             !-------------------------------
                             
                             IF ( ppm_kind == ppm_kind_double ) THEN
                                
                                ibuffer = ibuffer + 1 
                                ppm_sendbufferd(ibuffer) = xt(1,i)
                                ibuffer = ibuffer + 1 
                                ppm_sendbufferd(ibuffer) = xt(2,i)
                                ibuffer = ibuffer + 1 
                                ppm_sendbufferd(ibuffer) = xt(3,i)
                                
                             ELSE
                                
                                ibuffer = ibuffer + 1 
                                ppm_sendbuffers(ibuffer) = xt(1,i)
                                ibuffer = ibuffer + 1 
                                ppm_sendbuffers(ibuffer) = xt(2,i)
                                ibuffer = ibuffer + 1 
                                ppm_sendbuffers(ibuffer) = xt(3,i)
                                
                             END IF ! ppm_kind
                             
                          END IF ! xt
                        
                       END DO ! i =1, nghostplus
                       
                    ENDIF ! (2/3 dimension )
                    
                 ENDIF ! (only consider subs belonging to rank)
                 
              ENDDO ! j (loop over all subs in topo)
              
           ENDIF ! sendrank .GE. 0 (skipped negative ranks)
         
           !-------------------------------------------------
           ! Update the buffer pointer (i.e., the current iset
           ! or the number of particles we will send to the 
           ! k-th entry in the icommseq list
           !-------------------------------------------------
           
           ppm_psendbuffer(k+1) = iset + 1
           
        END DO ! k loop over all processors in commseq
        
        !----------------------------------------------------
        ! Store the current size of the buffer
        !----------------------------------------------------
        
        ppm_nsendbuffer = ibuffer
        
        !----------------------------------------------------
        ! Deallocate the memory for the lists
        !----------------------------------------------------
        
        iopt = ppm_param_dealloc
        CALL ppm_alloc(ilist1,ldu,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ghost_get', &
                'ilist1',__LINE__,info)
        END IF
        
        CALL ppm_alloc(ilist2,ldu,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ghost_get', &
                'ilist2',__LINE__,info)
        END IF
        
        CALL ppm_alloc(ighost,ldu,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ghost_get', &
                'ighost',__LINE__,info)
        END IF
        
        CALL ppm_alloc( xt,ldu,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ghost_get', &
                'xt',__LINE__,info)
        END IF
        
        CALL ppm_alloc( lghost,ldu,iopt,info)
        IF (info .NE. 0) THEN
           info = ppm_error_error
           CALL ppm_error(ppm_err_dealloc,'ppm_map_part_ghost_get', &
                'lghost',__LINE__,info)
        END IF
        
        !----------------------------------------------------
        ! Return 
        !----------------------------------------------------
        
9999    CONTINUE
        
        CALL substop('ppm_map_part_ghost_get',t0,info)
        
        RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_map_part_ghost_get_s
#elif  __KIND == __DOUBLE_PRECISION
    END SUBROUTINE ppm_map_part_ghost_get_d
#endif
    
