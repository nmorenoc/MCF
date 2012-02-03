#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_tree_cutdir
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine finds the n best directions of cut for
      !                 a given box.
      !
      !  Input        : xp(:,:)      (F) the positions of the particles
      !                 Npart        (I) the number of particles 
      !                 weights(3)   (F) weights for the three cost
      !                                  contributions: particles, mesh,
      !                                  geometry used to find the cut
      !                                  directions.
      !                 min_box(:,:) (F) the minimum coordinate of the 
      !                                  boxes
      !                 max_box(:,:) (F) the maximum coordinate of the 
      !                                  boxes
      !                 cutbox       (I) ID of the box to be cut.
      !                 ncut         (I) number of cut directions to be
      !                                  found
      !                 fixed(:)     (L) set to .TRUE. for dimensions which
      !                                  must not be cut
      !                 minboxsize(:)(F) minimum size of a box in all
      !                                  directions.
      !                 pcost(:)     (F) OPTIONAL argument of length
      !                                  Npart, specifying the
      !                                  computational cost of each
      !                                  particle.
      !
      !  Input/output :                                            
      !
      !  Output       : icut(:)      (I) directions of best cut. icut=i
      !                                  means: cutting plane is orthogonal
      !                                  to i-th coordinate axis. index:
      !                                  1..ncut. The directions are sorted
      !                                  with the most favorable first.
      !                 info         (I) return status
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_tree_cutdir.f,v $
      !  Revision 1.9  2006/09/05 08:01:27  pchatela
      !  Proper scaling for REAL comparisons
      !  Added module_alloc to ppm_decomp_boxsplit
      !
      !  Revision 1.8  2005/08/31 13:33:50  ivos
      !  bugfix: removed doubly-declared variables and unused arguments.
      !
      !  Revision 1.7  2005/08/31 11:24:31  ivos
      !  Further optimizations and bugfix in the maxcost computation.
      !
      !  Revision 1.6  2005/08/30 13:17:26  ivos
      !  Sharked the routines and unrolled all loops over ppm_dim.
      !
      !  Revision 1.5  2004/12/03 17:14:59  ivos
      !  Switched to the use of particle lists lhbx and lpdx.
      !
      !  Revision 1.4  2004/11/03 16:27:48  ivos
      !  Box is now centered at 0 and distances are scaled before computing the
      !  tensor of intertia. This improves the condition number of the
      !  matrix and prevents problems in the eigenvector computations.
      !
      !  Revision 1.3  2004/09/23 09:49:53  ivos
      !  bugfix: only directions which are long enough (.GT.2*minboxsize)
      !  are proposed for cutting.
      !
      !  Revision 1.2  2004/09/22 17:27:48  ivos
      !  bugfix: assignment of directions to icut was wrong.
      !
      !  Revision 1.1  2004/09/22 10:32:04  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_cutdir_s(xp,Npart,weights,min_box,max_box,   &
     &    cutbox,ncut,fixed,minboxsize,icut,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_cutdir_d(xp,Npart,weights,min_box,max_box,   &
     &    cutbox,ncut,fixed,minboxsize,icut,info,pcost)
#endif

      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_tree
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      USE ppm_module_util_eigen_2sym
      USE ppm_module_util_eigen_3sym
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
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
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_box,max_box
      REAL(MK), DIMENSION(3  ), INTENT(IN   ) :: minboxsize,weights
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      INTEGER                 , INTENT(IN   ) :: Npart,ncut,cutbox
      LOGICAL , DIMENSION(:  ), INTENT(IN   ) :: fixed
      INTEGER , DIMENSION(:  ), POINTER       :: icut
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)            :: len_box,EW
      REAL(MK), DIMENSION(ppm_dim,ppm_dim)    :: J,Jsum,EV
      REAL(MK), DIMENSION(ppm_dim)            :: cgeom,cmesh,cpart,ctotal
      REAL(MK), DIMENSION(ppm_dim)            :: shift
      INTEGER , DIMENSION(2)                  :: ldc
      INTEGER , DIMENSION(ppm_dim)            :: isort
      LOGICAL , DIMENSION(ppm_dim)            :: cutable
      REAL(MK)                                :: t0,csum,csuminv,lmyeps
      REAL(MK)                                :: a0,a1,a2,dis,x2,y2,z2
      REAL(MK)                                :: x,y,z,fscale
      INTEGER                                 :: i,k,ip,iopt
#ifdef USE_MPI
      INTEGER                                 :: MPTYPE
#endif
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree_cutdir',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (cutbox .LE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_cutdir',     &
     &          'cutbox must be > 0 !',__LINE__,info)
            GOTO 9999
         ENDIF 
         IF (SIZE(min_box,2) .LT. cutbox) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_cutdir',     &
     &          'size of min_box must be at least cutbox !',__LINE__,info)
            GOTO 9999
         ENDIF 
         IF (SIZE(max_box,2) .LT. cutbox) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_cutdir',     &
     &          'size of max_box must be at least cutbox !',__LINE__,info)
            GOTO 9999
         ENDIF 
         DO i=1,ppm_dim
            IF (min_box(i,cutbox) .GT. max_box(i,cutbox)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_tree_cutdir',     &
     &             'min_box must be <= max_box !',__LINE__,info)
               GOTO 9999
            ENDIF 
            IF (minboxsize(i) .LT. 0.0_MK) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_tree_cutdir',     &
     &             'minboxsize must be >= 0.0 !',__LINE__,info)
               GOTO 9999
            ENDIF 
         ENDDO
      ENDIF 

      !-------------------------------------------------------------------------
      !  If we have less than 1 direction to find, we are done
      !-------------------------------------------------------------------------
      IF (ncut .LT. 1) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree_cutdir',   &
     &            'No directions requested. Exiting.',info)
          ENDIF
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute the extension of the box
      !-------------------------------------------------------------------------
      IF (ppm_dim .GT. 2) THEN
          len_box(1) = max_box(1,cutbox)-min_box(1,cutbox)
          len_box(2) = max_box(2,cutbox)-min_box(2,cutbox)
          len_box(3) = max_box(3,cutbox)-min_box(3,cutbox)
      ELSE
          len_box(1) = max_box(1,cutbox)-min_box(1,cutbox)
          len_box(2) = max_box(2,cutbox)-min_box(2,cutbox)
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine cuttable directions of this box
      !-------------------------------------------------------------------------
      cutable = .FALSE.
      ip = 0
      DO i=1,ppm_dim
          IF ((.NOT.fixed(i)) .AND. (len_box(i)-(2.0_MK*minboxsize(i))   &
     &        .GT.lmyeps*len_box(i))) THEN
              ip = ip + 1
              cutable(i) = .TRUE.
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Exit if box has not enough cuttable directions
      !-------------------------------------------------------------------------
      IF (ip .LT. ncut) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree_cutdir',   &
     &            'Not enough cutable directions! Exiting.',info)
          ENDIF
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If this is an octtree in 3d or a quad-tree in 2d it is trivial
      !-------------------------------------------------------------------------
      IF (ncut .EQ. 3 .AND. ppm_dim .EQ. 3) THEN
          icut(1) = 1
          icut(2) = 2
          icut(3) = 3
          GOTO 9999
      ELSEIF (ncut .EQ. 2 .AND. ppm_dim .EQ. 2) THEN
          icut(1) = 1
          icut(2) = 2
          GOTO 9999
      ENDIF

#ifdef USE_MPI
      !-------------------------------------------------------------------------
      !  Determine MPI data type
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      MPTYPE = MPI_REAL
#elif __KIND == __DOUBLE_PRECISION 
      MPTYPE = MPI_DOUBLE_PRECISION
#endif
#endif
      !-------------------------------------------------------------------------
      !  Geometry: no mesh and no particles
      !-------------------------------------------------------------------------
      ctotal = 0.0_MK
      cgeom  = 0.0_MK
      IF (weights(3) .NE. 0.0_MK) THEN
          IF (ppm_dim .GT. 2) THEN
              csum = SUM(len_box(1:3))
          ELSE
              csum = SUM(len_box(1:2))
          ENDIF
          ! normalozation: total cost contribution of each part needs to
          ! sum up to 1 in order to get correct weighting.
          csuminv = 1.0_MK/csum
          DO i=1,ppm_dim
              cgeom(i)  = len_box(i)*csuminv
              ctotal(i) = weights(3)*cgeom(i)
          ENDDO
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Mesh costs
      !-------------------------------------------------------------------------
      cmesh = 0.0_MK
      IF (have_mesh .AND. weights(2) .NE. 0.0_MK) THEN
          IF (ppm_dim .GT. 2) THEN
              csum = SUM(REAL(Nm_box(1:3,cutbox),MK))
          ELSE
              csum = SUM(REAL(Nm_box(1:2,cutbox),MK))
          ENDIF
          ! normalozation: total cost contribution of each part needs to
          ! sum up to 1 in order to get correct weighting.
          csuminv = 1.0_MK/csum
          DO i=1,ppm_dim
              cmesh(i)  = REAL(Nm_box(i,cutbox),MK)*csuminv
              ctotal(i) = ctotal(i) + (weights(2)*cmesh(i))
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Particle costs
      !-------------------------------------------------------------------------
      cpart = 0.0_MK
      IF (have_particles .AND. weights(1) .NE. 0.0_MK) THEN
          !---------------------------------------------------------------------
          !  Shift box such that its center is at 0 and scale with max
          !  boxlength to get well-conditioned matrix
          !---------------------------------------------------------------------
          shift(1)   = 0.5_MK*(min_box(1,cutbox)+max_box(1,cutbox))
          shift(2)   = 0.5_MK*(min_box(2,cutbox)+max_box(2,cutbox))
          IF (ppm_dim .GT. 2) THEN
              shift(3)   = 0.5_MK*(min_box(3,cutbox)+max_box(3,cutbox))
              fscale = MAXVAL(len_box(1:3))
          ELSE
              fscale = MAXVAL(len_box(1:2))
          ENDIF
          fscale = 1.0_MK/fscale

          !---------------------------------------------------------------------
          !  Compute tensor of inertia of local particles
          !---------------------------------------------------------------------
          J = 0.0_MK
          IF (ppm_dim .EQ. 2) THEN
              DO k=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
                  ip = tree_lpdx(k)
                  x = (xp(1,ip) - shift(1))*fscale
                  y = (xp(2,ip) - shift(2))*fscale
                  x2 = x*x
                  y2 = y*y
                  IF (PRESENT(pcost)) THEN
                      J(1,1) = J(1,1) + (y2*pcost(ip))    ! Jx
                      J(2,2) = J(2,2) + (x2*pcost(ip))    ! Jy
                      J(1,2) = J(1,2) - (x*y*pcost(ip))   ! Cxy
                  ELSE
                      J(1,1) = J(1,1) + y2      ! Jx
                      J(2,2) = J(2,2) + x2      ! Jy
                      J(1,2) = J(1,2) - (x*y)   ! Cxy
                  ENDIF
              ENDDO
          ELSE
              DO k=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
                  ip = tree_lpdx(k)
                  x = (xp(1,ip) - shift(1))*fscale
                  y = (xp(2,ip) - shift(2))*fscale
                  z = (xp(3,ip) - shift(3))*fscale
                  x2 = x*x
                  y2 = y*y
                  z2 = z*z
                  IF (PRESENT(pcost)) THEN
                      J(1,1) = J(1,1) + ((y2+z2)*pcost(ip))    ! Jx
                      J(2,2) = J(2,2) + ((x2+z2)*pcost(ip))    ! Jy
                      J(3,3) = J(3,3) + ((x2+y2)*pcost(ip))    ! Jz
                      J(1,2) = J(1,2) - (x*y*pcost(ip))        ! Cxy
                      J(1,3) = J(1,3) - (x*z*pcost(ip))        ! Cxz
                      J(2,3) = J(2,3) - (y*z*pcost(ip))        ! Cyz
                  ELSE
                      J(1,1) = J(1,1) + (y2+z2)    ! Jx
                      J(2,2) = J(2,2) + (x2+z2)    ! Jy
                      J(3,3) = J(3,3) + (x2+y2)    ! Jz
                      J(1,2) = J(1,2) - (x*y)      ! Cxy
                      J(1,3) = J(1,3) - (x*z)      ! Cxz
                      J(2,3) = J(2,3) - (y*z)      ! Cyz
                  ENDIF
              ENDDO
          ENDIF

#ifdef USE_MPI
          !---------------------------------------------------------------------
          !  Allreduce of tensor of intertia
          !---------------------------------------------------------------------
          ip = 9
          IF (ppm_dim .EQ. 2) ip = 4
          CALL MPI_Allreduce(J,Jsum,ip,MPTYPE,MPI_SUM,ppm_comm,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_mpi_fail,'ppm_tree_cutdir',   &
     &            'MPI_Allreduce of inertia tensor failed!',__LINE__,info)
              GOTO 9999
          ENDIF 
          J = Jsum
#endif

          !---------------------------------------------------------------------
          !  Build full symmetric matrix for Eigendecomposition routine
          !---------------------------------------------------------------------
          J(2,1) = J(1,2)
          IF (ppm_dim .GT. 2) THEN
              J(3,1) = J(1,3)
              J(3,2) = J(2,3)
          ENDIF

          !---------------------------------------------------------------------
          !  Compute Eigenvectors of tensor of inertia
          !---------------------------------------------------------------------
          IF (ppm_dim .EQ. 3) THEN
              CALL ppm_util_eigen_3sym(J,EW,EV,info)
          ELSE
              CALL ppm_util_eigen_2sym(J,EW,EV,info)
          ENDIF
          IF (info.NE.0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_sub_failed,'ppm_tree_cutdir',   &
     &            'Eigendecomposition failed!',__LINE__,info)
              GOTO 9999
          ENDIF 

          !---------------------------------------------------------------------
          !  Flip Eigenvalues because we want the cost of a cut
          !  PERPENDICULAR to the EV directions and not parallel to them
          !---------------------------------------------------------------------
          csum        = EW(1)
          EW(1)       = EW(ppm_dim)
          EW(ppm_dim) = csum

          !---------------------------------------------------------------------
          !  Assign cost to coordinate directions
          !---------------------------------------------------------------------
          IF (ppm_dim .GT. 2) THEN
              csum = SUM(EW(1:3))
          ELSE
              csum = SUM(EW(1:2))
          ENDIF
          csuminv = 1.0_MK/csum
          DO i=1,ppm_dim
              DO k=1,ppm_dim
                  cpart(k) = cpart(k) + (EV(k,i)*EW(i)*csuminv)
              ENDDO
          ENDDO
          DO i=1,ppm_dim
              ctotal(i) = ctotal(i) + (weights(1)*cpart(i))
          ENDDO
      ENDIF     ! have_particles

      !-------------------------------------------------------------------------
      !  Sort directions by total cost. We do not need to normalize ctotal
      !  by SUM(weights), since we are only interested in the ORDER of
      !  directions and not the actual cost values.
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          isort(1:2) = (/1,2/)
          IF (ctotal(2) .GT. ctotal(1)) isort(1:2) = (/2,1/)
      ELSE
          isort(1:3) = (/1,2,3/)
          IF (ctotal(2) .GT. ctotal(1)) isort(1:3) = (/2,1,3/)
          IF (ctotal(3).GT.ctotal(1).AND.ctotal(3).GT.ctotal(2)) THEN
              isort(3) = isort(2)
              isort(2) = isort(1)
              isort(1) = 3
          ELSEIF (ctotal(3).GT.ctotal(isort(2))) THEN
              isort(3) = isort(2)
              isort(2) = 3
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Assign ncut best non-fixed directions in descending order to icut
      !-------------------------------------------------------------------------
      ip = 0
      i  = 1
      DO WHILE ((ip .LT. ncut) .AND. (i .LE. ppm_dim))
          k = isort(i)
          IF (cutable(k)) THEN
              ip = ip + 1
              icut(ip) = k
          ENDIF
          i = i + 1
      ENDDO

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_tree_cutdir',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_cutdir_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_cutdir_d
#endif
