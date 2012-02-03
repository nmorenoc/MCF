#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_template_comp_pp_ring
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Template for the user subroutine which computes
      !                 direct particle-particle interactions using the
      !                 ring topology. How to use this template:
      !                    (1) copy this file to your application code dir
      !                    (2) rename the file and this routine (make sure
      !                        to also replace all the occurrences of the
      !                        subroutine name in all WRITE statements).
      !                    (3) change address in header and add your CVS
      !                        Log (if needed).
      !                    (4) declare and add additional arguments
      !                    (5) add USE and INCLUDE statements for your
      !                        modules and header files
      !                    (6) define the precision of floating point data
      !                        (maybe provide different versions of this
      !                        subroutine for different precisions)
      !                    (7) add particle-particle interactions where
      !                        needed (search for USER CODE)
      !                    (8) Optional: Make argument checking conditional
      !                        on the debug level of your application
      !                    (9) Optional: Make this routine use your own
      !                        error logging, message writing, allocation,
      !                        subroutine start, subroutine stop and
      !                        timing routines.
      !
      !  Input        : xp(:,:)      (F) : particle co-ordinates [Group 1]
      !                 vp(:,:)      (F) : particle data used for interaction
      !                                    (e.g. vorticity, strength, ...)
      !                                    [Group 1]
      !                 ...          (.) : USER: add more input arguments for
      !                                    group 1 here if needed.
      !                                    needed
      !                 Npart        (I) : number of particles on this proc.
      !                                    [Group 1]
      !                 xp2(:,:)     (F) : particle co-ordinates [Group 2]
      !                 vp2(:,:)     (F) : particle data used for interaction
      !                                    (e.g. vorticity, strength, ...)
      !                                    [Group 2]
      !                 ...          (.) : USER: add more input arguments for
      !                                    group 1 here if needed.
      !                 Lpart        (I) : number of particles [Group 2]
      !                 lsymm        (L) : Whether to use symmetry or not:
      !                                     .FALSE. pp interaction w/o symmetry
      !                                     .TRUE.  pp interaction w/  symmetry
      !                 mode         (I) : Whether the two groups are the
      !                                    same or not:
      !                                     0 not the same group
      !                                     1 the same group
      !
      !  Input/output : params(:)    (F) : user defined parameters and/or
      !                                    output from the interaction (i.e.
      !                                    potential energy)
      !                 fp(:,:)      (F) : Change of particle data due
      !                                    to interaction [Group 1]
      !                 fp2(:,:)     (F) : Change of particle data due
      !                                    to interaction [Group 2]
      !                 ...          (.) : USER: add more output arguments for
      !                                    if needed.
      !  Output       : info         (I) : return status. 0 if no error.
      !
      !  Routines     : 
      !
      !  Remarks      : Search for USER in this file and fill in
      !                 the particle-particle interactions in all 4 places.
      !                 Use the fp and fp2 array to return the result.
      !
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Revision: 1.2 $
      !  Revision 1.2  2004/04/23 17:25:33  oingo
      !  Revision 1.1  2004/04/22 08:29:06  oingo
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_template_comp_pp_ring(xp, vp, fp ,Npart,        &
     &                                     xp2,vp2,fp2,Lpart,        &
     &                                     lsymm,params,mode,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      ! USER: USE statements here...
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      INCLUDE 'ppm_param.h'
      !-------------------------------------------------------------------------
      !  USER: Define precision
      !-------------------------------------------------------------------------
!     INTEGER, PARAMETER :: MK = ppm_kind_single
      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: vp
      ! USER: add more arguments here if needed
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: fp
      INTEGER                 , INTENT(IN   ) :: Npart
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp2
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: vp2
      ! USER: add more arguments here if needed
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: fp2
      INTEGER                 , INTENT(IN   ) :: Lpart
      LOGICAL                 , INTENT(IN   ) :: lsymm
      REAL(MK), DIMENSION(:),   INTENT(INOUT) :: params
      INTEGER                 , INTENT(IN   ) :: mode
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      ! counters
      INTEGER                                    :: i,j
      ! coordinate differences
      REAL(MK)                                   :: dx,dy,dz
      ! square of inter particle distance
      REAL(MK)                                   :: dij

      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      info = 0

      !-------------------------------------------------------------------------
      !  Check arguments
      !  USER: You might want to make this conditional on the debug level
      !  of your application...
      !-------------------------------------------------------------------------
      IF (Npart .LE. 0) THEN
         WRITE(*,'(2A)') '(ppm_template_comp_pp_ring): ',  &
     &        'Npart must be >0'
         info = -1
         GOTO 9999
      ENDIF
      IF (Lpart .LE. 0) THEN
         WRITE(*,'(2A)') '(ppm_template_comp_pp_ring): ',  &
     &        'Lpart must be >0'
         info = -1
         GOTO 9999
      ENDIF
      IF ((mode .NE. 0) .AND. (mode .NE. 1)) THEN
         WRITE(*,'(2A)') '(ppm_template_comp_pp_ring): ',  &
     &        'MODE must be either 0 or 1'
         info = -1
         GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Check if we are computing the selfinteraction (i.e. Group1 .EQ. Group2)
      !-------------------------------------------------------------------------
      IF (mode .EQ. 1) THEN
         IF (lsymm) THEN
            DO i = 1,Npart
               DO j = i+1,Lpart
                  !-----------------------------------------------------------
                  !  USER CODE HERE (with symmetry)
                  !  dx = xp(1,i) - xp2(1,j)
                  !  dy = xp(2,i) - xp2(2,j)
                  !  dz = xp(3,i) - xp2(3,j)
                  !  dij = (dx*dx)+(dy*dy)+(dz*dz)
                  !  ...
                  !-----------------------------------------------------------
               ENDDO
            ENDDO
         ELSE
            DO i = 1,Npart
               DO j = 1,Lpart
                  IF (i .EQ. j) CYCLE
                  !-----------------------------------------------------------
                  !  USER CODE HERE (without symmetry)
                  !-----------------------------------------------------------
               ENDDO
            ENDDO
         ENDIF
      !-------------------------------------------------------------------------
      !  Here we compute the interaction between two different groups
      !-------------------------------------------------------------------------
      ELSE
         IF (lsymm) THEN
            DO i = 1,Npart
               DO j = 1,Lpart
                  !-----------------------------------------------------------
                  !  USER CODE HERE (with symmetry)
                  !-----------------------------------------------------------
               ENDDO
            ENDDO
         ELSE
            DO i = 1,Npart
               DO j = 1,Lpart
                  !-----------------------------------------------------------
                  !  USER CODE HERE (without symmetry)
                  !-----------------------------------------------------------
               ENDDO
            ENDDO
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_template_comp_pp_ring
