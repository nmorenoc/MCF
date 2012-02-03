#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_tree_cutpos
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine finds the best cutting positions for
      !                 the given cut directions.
      !
      !  Input        : xp(:,:)      (F) the positions of the particles
      !                 Npart        (I) the number of particles 
      !                 weights(3)   (F) weights for the three cost
      !                                  contributions: particles, mesh,
      !                                  geometry for finding the cut
      !                                  positions.
      !                 min_box(:,:) (F) the minimum coordinate of the 
      !                                  boxes
      !                 max_box(:,:) (F) the maximum coordinate of the 
      !                                  boxes
      !                 cutbox       (I) ID of box to be cut
      !                 ncut         (I) number of cut directions 
      !                 minboxsize(:)(F) minimum box size required in all
      !                                  spatial directions
      !                 icut(:)      (I) cut directions
      !                 pcost(:)     (F) OPTIONAL argument of length
      !                                  Npart, specifying the
      !                                  computational cost of each
      !                                  particle.
      !
      !  Input/output :                                            
      !
      !  Output       : cpos(:)      (F) positions of best cuts. index:
      !                                  1..ncut. 
      !                 info         (I) return status
      !
      !  Remarks      : The cost of a particle is counted as pcost (if
      !                 given) or 1. For meshes, the cost is 1 per mesh
      !                 point. For the geometry part, the cost of a
      !                 box is given by its volume. Use weights(3) to
      !                 adjust this if needed.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_tree_cutpos.f,v $
      !  Revision 1.11  2005/08/31 13:33:50  ivos
      !  bugfix: removed doubly-declared variables and unused arguments.
      !
      !  Revision 1.10  2005/08/31 11:24:32  ivos
      !  Further optimizations and bugfix in the maxcost computation.
      !
      !  Revision 1.9  2005/08/30 13:17:27  ivos
      !  Sharked the routines and unrolled all loops over ppm_dim.
      !
      !  Revision 1.8  2004/12/03 17:15:00  ivos
      !  Switched to the use of particle lists lhbx and lpdx.
      !
      !  Revision 1.7  2004/12/02 16:31:54  ivos
      !  Optimized some loops.
      !
      !  Revision 1.6  2004/12/02 10:00:38  ivos
      !  bugfix: total cost per box was counted ncut times.
      !
      !  Revision 1.5  2004/09/30 10:21:43  ivos
      !  bugfix: in the case without particles, cp was 0 and partpos thus
      !  NAN. Fixed it by adding IF statement around partpos calculation.
      !
      !  Revision 1.4  2004/09/24 08:01:27  ivos
      !  bugfix: index error fixed.
      !
      !  Revision 1.3  2004/09/23 09:49:06  ivos
      !  fixed indentation.
      !
      !  Revision 1.2  2004/09/22 17:27:08  ivos
      !  bugfix: added USE ppm_module_write
      !
      !  Revision 1.1  2004/09/22 10:32:05  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_cutpos_s(xp,Npart,weights,min_box,max_box,   &
     &    cutbox,ncut,minboxsize,icut,cpos,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_cutpos_d(xp,Npart,weights,min_box,max_box,   &
     &    cutbox,ncut,minboxsize,icut,cpos,info,pcost)
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
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: icut
      REAL(MK), DIMENSION(:  ), POINTER       :: cpos
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)            :: len_box
      REAL(MK), DIMENSION(ncut+1)             :: pc,pcsum
      INTEGER , DIMENSION(2)                  :: ldc
      REAL(MK)                                :: t0,dm,meshtotal,geomtotal
      REAL(MK)                                :: pmass,mmass,gmass,tmass
      REAL(MK)                                :: partpos,midpos
      INTEGER                                 :: i,j,ip,cutdir,ncp1,iopt
      INTEGER                                 :: MPTYPE
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree_cutpos',t0,info)

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         IF (cutbox .LE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_cutpos',     &
     &          'cutbox must be > 0 !',__LINE__,info)
            GOTO 9999
         ENDIF 
         IF (SIZE(min_box,2) .LT. cutbox) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_cutpos',     &
     &          'size of min_box must be at least cutbox !',__LINE__,info)
            GOTO 9999
         ENDIF 
         IF (SIZE(max_box,2) .LT. cutbox) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_cutpos',     &
     &          'size of max_box must be at least cutbox !',__LINE__,info)
            GOTO 9999
         ENDIF 
         DO i=1,ppm_dim
            IF (min_box(i,cutbox) .GT. max_box(i,cutbox)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_tree_cutpos',     &
     &             'min_box must be <= max_box !',__LINE__,info)
               GOTO 9999
            ENDIF 
         ENDDO
      ENDIF 

      !-------------------------------------------------------------------------
      !  If we have less than 1 direction to cut, we are done
      !-------------------------------------------------------------------------
      IF (ncut .LT. 1) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree_cutpos',   &
     &            'No cut directions present. Exiting.',info)
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

      pc = 0.0_MK
      ncp1 = ncut+1
      !-------------------------------------------------------------------------
      !  Project particle cost onto cut directions
      !-------------------------------------------------------------------------
      IF (have_particles .AND. weights(1) .NE. 0.0_MK) THEN
#ifdef __VECTOR
          DO i=1,ncut
              cutdir = icut(i)
              DO j=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
                  ip = tree_lpdx(j)
                  dm = 1.0_MK
                  IF (PRESENT(pcost)) dm = pcost(ip)
                  pc(i) = pc(i) + (xp(cutdir,ip)*dm)
                  pc(ncp1) = pc(ncp1) + dm
              ENDDO
          ENDDO
          ! replace this division by sth more clever! Count correctly in
          ! the first place
          pc(ncp1) = pc(ncp1)/REAL(ncut,MK)
#else
          DO j=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
              ip = tree_lpdx(j)
              dm = 1.0_MK
              IF (PRESENT(pcost)) dm = pcost(ip)
              DO i=1,ncut
                  cutdir = icut(i)
                  pc(i) = pc(i) + (xp(cutdir,ip)*dm)
              ENDDO
              pc(ncp1) = pc(ncp1) + dm
          ENDDO
#endif

#ifdef USE_MPI
          !---------------------------------------------------------------------
          !  Allreduce of projected particle sums
          !---------------------------------------------------------------------
          CALL MPI_Allreduce(pc,pcsum,ncp1,MPTYPE,MPI_SUM,ppm_comm,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_mpi_fail,'ppm_tree_cutpos',   &
     &            'MPI_Allreduce of projected particles',__LINE__,info)
              GOTO 9999
          ENDIF 
          pc = pcsum
#endif
      ENDIF   ! have_particles

      !-------------------------------------------------------------------------
      !  Total weight of mesh and geometry. 
      !  Convert to REAL first to avoid integer overflow.
      !-------------------------------------------------------------------------
      meshtotal = 0.0_MK
      IF (ppm_dim .EQ. 2) THEN
          IF (have_mesh .AND. weights(2) .NE. 0) THEN
              meshtotal = REAL(Nm_box(1,cutbox),MK)*REAL(Nm_box(2,cutbox),MK)
          ENDIF
          geomtotal = len_box(1)*len_box(2)
      ELSE
          IF (have_mesh .AND. weights(2) .NE. 0) THEN
              meshtotal = REAL(Nm_box(1,cutbox),MK)*    &
     &            REAL(Nm_box(2,cutbox),MK)*REAL(Nm_box(3,cutbox),MK)
          ENDIF
          geomtotal = len_box(1)*len_box(2)*len_box(3)
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute weighted masses of particles, mesh and geometry
      !-------------------------------------------------------------------------
      pmass = pc(ncp1) *weights(1)
      mmass = meshtotal*weights(2)
      gmass = geomtotal*weights(3)
      tmass = pmass+mmass+gmass

      IF (tmass .EQ. 0.0_MK) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_argument,'ppm_tree_cutpos',     &
     &        'Total cost is 0! Are all weights 0?',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  The optimal cut position is in the weighted center of mass
      !-------------------------------------------------------------------------
      DO i=1,ncut
          cutdir  = icut(i)
          IF (have_particles .AND. weights(1) .NE. 0.0_MK) THEN
              partpos = pc(i)/pc(ncp1)
          ELSE
              partpos = 0.0_MK
          ENDIF
          midpos  = min_box(cutdir,cutbox) + (0.5_MK*len_box(cutdir))
          cpos(i) = partpos*pmass + midpos*(mmass+gmass)
          cpos(i) = cpos(i)/tmass

          !---------------------------------------------------------------------
          !  Enforce that minboxsize is respected.
          !---------------------------------------------------------------------
          IF (cpos(i)-min_box(cutdir,cutbox) .LT. minboxsize(cutdir)) THEN
              cpos(i) = min_box(cutdir,cutbox)+minboxsize(cutdir)
          ENDIF
          IF (max_box(cutdir,cutbox)-cpos(i) .LT. minboxsize(cutdir)) THEN
              cpos(i) = max_box(cutdir,cutbox)-minboxsize(cutdir)
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_tree_cutpos',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_cutpos_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_cutpos_d
#endif
