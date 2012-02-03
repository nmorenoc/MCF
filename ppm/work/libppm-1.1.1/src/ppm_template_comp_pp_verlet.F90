#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :            ppm_template_comp_pp_verlet
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Template for the user subroutine which computes
      !                 direct particle-particle interactions using Verlet
      !                 lists. How to use this template:
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
      !                        needed (search for USER CODE). Make sure to
      !                        replace dimens with the variable which
      !                        actually contains the problem
      !                        dimensionality in your code.
      !                    (8) Optional: Make argument checking conditional
      !                        on the debug level of your application
      !                    (9) Optional: Make this routine use your own
      !                        error logging, message writing, allocation,
      !                        subroutine start, subroutine stop and
      !                        timing routines.
      !                    (10)Optional: Move the variable declarations for
      !                        vlist and nvlist to your client program and
      !                        add them as arguments to this routine. This
      !                        makes it possible to use several sets of
      !                        lists. 
      !                    (11)Optional: Remove the imode=-1 and imode=1
      !                        cases from this routine and build/destroy
      !                        the lists directly in your client program by
      !                        calling ppm_neighlist_vlist/DEALLOCATE. You
      !                        can then remove imode from the argument list
      !                        of this routine.
      !
      !  Input        : xp(:,:)    (F) particle co-ordinates
      !                 pdata(:,:) (F) particle data used for interaction
      !                                (e.g. vorticity, strength, ...)
      !                 ...        (.) USER: add more input arguments as
      !                                      needed
      !                 Np         (I) number of particles on this proc.
      !                 cutoff     (F) cutoff radius for PP interactions
      !                 skin       (F) skin thikness for the Verlet list
      !                 lsymm      (L) use symmetry for PP interactions?
      !                                   .TRUE. : Yes
      !                                   .FALSE.: No
      !                 imode      (I) Mode of action. Any of the folowing:
      !                                 0  PP interactions 
      !                                 1  build Verlet lists and return
      !                                 -1 destroy Verlet lists and return
      !
      !  Input/output : 
      !
      !  Output       : dpd(:,:)   (F) Change of particle data (pdata) due to 
      !                                interaction.
      !                 info       (I) return status. =0 if no error.
      !
      !  Routines     : ppm_neighlist_vlist
      !
      !  Remarks      : If particles have moved by more than the skin
      !                 thikness between calls or lsymm has
      !                 changed, always call with imode=1.
      !
      !                 seach for USER in this file and fill in
      !                 the particle-particle interactions in both places.
      !                 Use the dpd array to return the result. It is your
      !                 responsibility to properly allocate and initialize 
      !                 this array BEFORE calling this routine.
      !
      !                 After finishing all time steps, always call this
      !                 routine with imode=-1 to free the memory of the
      !                 cell lists. Failure to do so will result in a
      !                 memory leak in your program.
      !
      !                 If the CYCLE command is a problem on your hardware
      !                 (e.g. prevents vectorization), split the DO-loop in
      !                 two parts as:
      !                                 DO jpart = istart,ipart-1
      !                                 ...
      !                                 ENDDO
      !                                 DO jpart = ipart+1,iend
      !                                 ...
      !                                 ENDDO
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Revision: 1.3 $
      !  Revision 1.2  2004/02/24 11:35:54  ivos
      !  Revision 1.1  2004/01/26 17:24:35  ivos
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_template_comp_pp_verlet(xp,pdata,Np,cutoff,skin,   &
     &               lsymm,imode,dpd,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      ! USER: USE statements here...
      USE ppm_module_neighlist
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
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: pdata
      ! USER: add more input arguments as needed
      INTEGER                 , INTENT(IN   ) :: Np
      REAL(MK)                , INTENT(IN   ) :: cutoff
      REAL(MK)                , INTENT(IN   ) :: skin
      LOGICAL                 , INTENT(IN   ) :: lsymm
      INTEGER                 , INTENT(IN   ) :: imode
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: dpd
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! Verlet lists
      ! USER: if allocation of the following two varaibles fails due to stack 
      ! size limitations, try putting them in a module and add a USE
      ! statement for it above. In order to use several Verlet lists, move
      ! these declarations to the client program and add vlist and nvlist
      ! to the argument list of this routine.
      INTEGER, DIMENSION(:,:), POINTER, SAVE     :: vlist
      INTEGER, DIMENSION(  :), POINTER, SAVE     :: nvlist
      ! counters
      INTEGER                                    :: jpart,ip,jp
      ! coordinate differences
      REAL(MK)                                   :: dx,dy,dz
      ! square of inter particle distance
      REAL(MK)                                   :: dij
      ! cutoff squared
      REAL(MK)                                   :: cut2
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      info = 0
      cut2 = cutoff*cutoff

      !-------------------------------------------------------------------------
      !  Check arguments.
      !  USER: You might want to make this conditional on the debug level of
      !  your application...
      !-------------------------------------------------------------------------
      IF (cutoff .LE. 0.0_MK) THEN
          WRITE(*,'(2A)') '(ppm_template_comp_pp_verlet): ',  &
     &        'cutoff must be >0'
          info = -1
          GOTO 9999
      ENDIF
      IF (skin .LT. 0.0_MK) THEN
          WRITE(*,'(2A)') '(ppm_template_comp_pp_verlet): ',  &
     &        'skin must be >0'
          info = -1
          GOTO 9999
      ENDIF
      IF (Np .LE. 0) THEN
          WRITE(*,'(2A)') '(ppm_template_comp_pp_verlet):',  &
     &        'Np must be >0'
          info = -1
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If imode = -1, destroy Verlet lists and return.
      !  USER: This could also be done by your client directly. The code is
      !  just here as a template.
      !-------------------------------------------------------------------------
      IF (imode .EQ. -1) THEN
          DEALLOCATE(vlist,nvlist,STAT=info)
          IF (info .NE. 0) THEN
              WRITE(*,'(2A,I8)') '(ppm_template_comp_pp_verlet): ',  &
     &            'DEALLOCATE failed on line ',__LINE__
              info = -1
          ENDIF  
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Build new Verlet lists if needed.
      !  USER: This could also be done by your client program calling
      !  ppm_neighlist_vlist directly. The code is just here as a template.
      !-------------------------------------------------------------------------
      IF (imode .EQ. 1) THEN
          !---------------------------------------------------------------------
          !  Generate Verlet lists
          !---------------------------------------------------------------------
          CALL ppm_neighlist_vlist(xp,Np,cutoff,skin,lsymm,vlist,nvlist,info)
          IF (info .NE. 0) THEN
              WRITE(*,'(2A,I8)') '(ppm_template_comp_pp_verlet): ',&
     &            'Building Verlet lists failed on line ',__LINE__
              info = -1
          ENDIF
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  PARTICLE-PARTICLE INTERACTIONS using symmetry
      !-------------------------------------------------------------------------
      IF (lsymm) THEN
          DO ip=1,Np
              DO jpart=1,nvlist(ip)
                  jp = vlist(jpart,ip)
                  !-------------------------------------------------------------
                  !  Calculate the square of the distance between the two 
                  !  particles. It will always be .LE. (cutoff+skin)**2 by 
                  !  construction of the Verlet list.
                  !  COMMENT THIS IF YOU DO NOT NEED THE DISTANCE!
                  !-------------------------------------------------------------
                  dx  = xp(1,ip) - xp(1,jp)
                  dy  = xp(2,ip) - xp(2,jp)
                  IF (dimens .GT. 2) THEN
                      dz  = xp(3,ip) - xp(3,jp)
                      dij = (dx*dx) + (dy*dy) + (dz*dz)
                  ELSE
                      dz = 0.0_MK
                      dij = (dx*dx) + (dy*dy)
                  ENDIF
                  ! skip this interaction if the particles are further
                  ! apart than the given cutoff
                  IF (dij .GT. cut2) CYCLE

                  !-------------------------------------------------------------
                  ! Particle ip interacts with particle jp here... and
                  ! vice versa to use symmetry.
                  !-------------------------------------------------------------
                  !
                  ! USER CODE HERE:
                  ! dpd(:,ip) = f(pdata(:,ip),pdata(:,jp),dx,dy,dz)
                  ! dpd(:,jp) = -dpd(:,ip)
              ENDDO
          ENDDO

      !-------------------------------------------------------------------------
      !  PARTICLE-PARTICLE INTERACTIONS not using symmetry
      !  This will not vectorize and has therefore been moved into a
      !  separate loop !!
      !-------------------------------------------------------------------------
      ELSE
          DO ip=1,Np
              DO jpart=1,nvlist(ip)
                  jp = vlist(jpart,ip)
                  !-------------------------------------------------------------
                  !  Calculate the square of the distance between the two 
                  !  particles. It will always be .LE. (cutoff+skin)**2 by 
                  !  construction of the Verlet list.
                  !  COMMENT THIS IF YOU DO NOT NEED THE DISTANCE!
                  !-------------------------------------------------------------
                  dx  = xp(1,ip) - xp(1,jp)
                  dy  = xp(2,ip) - xp(2,jp)
                  IF (dimens .GT. 2) THEN
                      dz  = xp(3,ip) - xp(3,jp)
                      dij = (dx*dx) + (dy*dy) + (dz*dz)
                  ELSE
                      dz = 0.0_MK
                      dij = (dx*dx) + (dy*dy)
                  ENDIF
                  ! skip this interaction if the particles are further
                  ! apart than the given cutoff
                  IF (dij .GT. cut2) CYCLE

                  !-------------------------------------------------------------
                  ! Particle ip interacts with particle jp.
                  !-------------------------------------------------------------
                  !
                  ! USER CODE HERE:
                  ! dpd(:,ip) = f(pdata(:,ip),pdata(:,jp),dx,dy,dz
              ENDDO
          ENDDO
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_template_comp_pp_verlet
