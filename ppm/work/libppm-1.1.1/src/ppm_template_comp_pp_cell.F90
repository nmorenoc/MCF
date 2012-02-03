#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_template_comp_pp_cell
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Template for the user subroutine which computes
      !                 direct particle-particle interactions using cell
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
      !                        replace all occurrences of dimens with the
      !                        variable which actually contains the problem
      !                        dimensionality in your code.
      !                    (8) Optional: Make argument checking conditional
      !                        on the debug level of your application
      !                    (9) Optional: Make this routine use your own
      !                        error logging, message writing, allocation,
      !                        subroutine start, subroutine stop and
      !                        timing routines.
      !                    (10)Optional: Move the variable declarations for
      !                        clist and Nm to your client program and
      !                        add them as arguments to this routine. This
      !                        makes it possible to use several sets of
      !                        lists. 
      !                    (11)Optional: Remove the imode=-1 and imode=1
      !                        cases from this routine and build/destroy
      !                        the lists directly in your client program by
      !                        calling ppm_neighlist_clist/..destroy. You
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
      !                 lsymm      (L) Do symmetric PP interactions?
      !                                   .TRUE.  : Yes
      !                                   .FALSE. : No
      !                 imode      (I) Mode of action. Any of the folowing:
      !                                 0  PP interactions
      !                                 1  build cell lists and return
      !                                 -1 destroy cell lists and return
      !                 nsublist   (I) number of subdomains on the local
      !                                processor
      !
      !  Input/output : params(:)  (F) user defined parameters and/or
      !                                output from the interaction (i.e.
      !                                potential energy)
      !
      !  Output       : dpd(:,:)   (F) Change of particle data (pdata) due to 
      !                                interaction.
      !                 info       (I) return status. =0 if no error.
      !
      !  Routines     : ppm_neighlist_clist
      !                 ppm_clist_destroy
      !                 ppm_neighlist_MkNeighIdx
      !
      !  Remarks      : If particles have moved between calls or lsymm has
      !                 changed, always call with imode=1.
      !
      !                 seach for USER in this file and fill in
      !                 the particle-particle interactions in all 4 places.
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
      !  $Revision: 1.5 $
      !  Revision 1.4  2004/02/24 11:35:53  ivos
      !  Revision 1.3  2004/02/19 14:01:39  gonnetp
      !  Revision 1.2  2004/02/04 17:20:48  ivos
      !  Revision 1.1  2004/01/26 17:24:35  ivos
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_template_comp_pp_cell(xp,pdata,Np,cutoff,lsymm,   &
     &               imode,nsublist,dpd,params,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      ! USER: USE modules here
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
      LOGICAL                 , INTENT(IN   ) :: lsymm
      INTEGER                 , INTENT(IN   ) :: imode,nsublist
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: dpd
      REAL(MK), DIMENSION(:),   INTENT(INOUT) :: params
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! counters
      INTEGER                                    :: i,idom,ibox,jbox
      INTEGER                                    :: ipart,jpart,ip,jp,isize
      INTEGER                                    :: cbox,iinter,j,k
      ! coordinate differences
      REAL(MK)                                   :: dx,dy,dz
      ! square of inter particle distance
      REAL(MK)                                   :: dij
      ! cutoff squared
      REAL(MK), SAVE                             :: cut2
      ! start and end particle in a box
      INTEGER                                    :: istart,iend,jstart,jend
      ! box size for cell list
      REAL(MK), DIMENSION(3)                     :: bsize
      ! cell list
      ! USER: if allocation of the following varaible fails due to stack 
      ! size limitations, try putting it in a module and add a USE
      ! statement for it above. In order to use several Cell lists, move
      ! these declarations to the client program and add clist and Nm
      ! to the argument list of this routine.
      TYPE(ppm_type_ptr_to_clist), DIMENSION(:), POINTER, SAVE :: clist
      ! number of cells in all directions
      INTEGER, DIMENSION(:,:), POINTER, SAVE     :: Nm
      ! cell neighbor lists
      INTEGER, DIMENSION(:,:), POINTER, SAVE     :: inp,jnp
      ! number of interactions for each cell
      INTEGER, SAVE                              :: nnp
      ! cell offsets for box index
      INTEGER                                    :: n1,n2,nz
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      info = 0

      !-------------------------------------------------------------------------
      !  Check arguments.
      !  USER: you may want to make this conditional on the debug level of
      !  your application...
      !-------------------------------------------------------------------------
      IF (cutoff .LE. 0.0_MK) THEN
          WRITE(*,'(2A)') '(ppm_template_comp_pp_cell): ',  &
     &        'cutoff must be >0'
          info = -1
          GOTO 9999
      ENDIF
      IF (Np .LE. 0) THEN
          WRITE(*,'(2A)') '(ppm_template_comp_pp_cell): ',  &
     &        'Np must be >0'
          info = -1
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If imode = -1, destroy cell lists and return
      !  USER: This could also be done by your client directly. The code is
      !  just here as a template.
      !-------------------------------------------------------------------------
      IF (imode .EQ. -1) THEN
          CALL ppm_clist_destroy(clist,info)
          DEALLOCATE(inp,jnp,Nm,STAT=info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Build new cell lists if needed
      !  USER: This could also be done by your client directly. The code is
      !  just here as a template.
      !-------------------------------------------------------------------------
      IF (imode .EQ. 1) THEN
          !---------------------------------------------------------------------
          !  Destroy old cell list (if there already was one)
          !---------------------------------------------------------------------
          CALL ppm_clist_destroy(clist,info)

          !---------------------------------------------------------------------
          !  The size of the cells in each spatial direction should at last
          !  be the size of the cutoff.
          !---------------------------------------------------------------------
          bsize(1:3) = cutoff
          cut2 = cutoff*cutoff

          !---------------------------------------------------------------------
          !  Generate cell lists
          !---------------------------------------------------------------------
          CALL ppm_neighlist_clist(xp,Np,bsize,lsymm,clist,Nm,info)
          IF (info .NE. 0) THEN
              WRITE(*,'(2A)') '(ppm_template_comp_pp_cell): ', &
     &            'Building cell lists failed '
              info = -1
              GOTO 9999
          ENDIF

          !---------------------------------------------------------------------
          !  Generate cell neighbor lists 
          !---------------------------------------------------------------------
          CALL ppm_neighlist_MkNeighIdx(lsymm,inp,jnp,nnp,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  PARTICLE-PARTICLE INTERACTIONS using symmetry
      !  This will not vectorize and has therefore been moved into a
      !  separate loop !!
      !-------------------------------------------------------------------------
      IF (lsymm) THEN
          DO idom=1,nsublist
              n1  = Nm(1,idom)
              n2  = Nm(1,idom)*Nm(2,idom)
              nz  = Nm(3,idom)
              IF (dimens .EQ. 2) THEN
                  n2 = 0
                  nz = 2
              ENDIF 
              ! loop over all REAL cells (the -2 in the end does this)
              DO k=0,nz-2
                  DO j=0,Nm(2,idom)-2
                      DO i=0,Nm(1,idom)-2
                          ! index of the center box
                          cbox = i + 1 + n1*j + n2*k
                          ! loop over all box-box interactions
                          DO iinter=1,nnp
                              ! determine box indices for this interaction
                              ibox = cbox+(inp(1,iinter)+n1*inp(2,iinter)+ &
     &                               n2*inp(3,iinter))
                              jbox = cbox+(jnp(1,iinter)+n1*jnp(2,iinter)+ &
     &                               n2*jnp(3,iinter))
                              !-------------------------------------------------
                              !  Read indices and check if cell is empty
                              !-------------------------------------------------
                              istart = clist(idom)%lhbx(ibox)
                              iend   = clist(idom)%lhbx(ibox+1)-1
                              IF (iend .LT. istart) CYCLE
                              !-------------------------------------------------
                              !  Within the box itself use symmetry and avoid 
                              !  adding the particle itself to its own list
                              !-------------------------------------------------
                              IF (ibox .EQ. jbox) THEN
                                  DO ipart=istart,iend-1
                                      ip = clist(idom)%lpdx(ipart)
                                      DO jpart=(ipart+1),iend
                                          jp = clist(idom)%lpdx(jpart)
                                          dx  = xp(1,ip) - xp(1,jp)
                                          dy  = xp(2,ip) - xp(2,jp)
                                          IF (dimens .GT. 2) THEN
                                              dz  = xp(3,ip) - xp(3,jp)
                                              dij = (dx*dx)+(dy*dy)+(dz*dz)
                                          ELSE
                                              dz = 0.0_MK
                                              dij = (dx*dx)+(dy*dy)
                                          ENDIF
                                          IF (dij .GT. cut2) CYCLE
                                          !-------------------------------------
                                          ! Particle ip interacts with
                                          ! particle jp here... and
                                          ! vice versa to use symmetry.
                                          !-------------------------------------
                                          !
                                          ! USER CODE HERE:
                                          ! dpd(:,ip) = f(pdata(:,ip),
                                          !               pdata(:,jp),dij)
                                          ! dpd(:,jp) = -dpd(:,ip)
                                      ENDDO
                                  ENDDO
                              !-------------------------------------------------
                              !  For the other boxes check all particles
                              !-------------------------------------------------
                              ELSE
                                  ! get pointers to first and last particle 
                                  jstart = clist(idom)%lhbx(jbox)
                                  jend   = clist(idom)%lhbx(jbox+1)-1
                                  ! skip this iinter if empty
                                  IF (jend .LT. jstart) CYCLE
                                  ! loop over all particles inside this cell
                                  DO ipart=istart,iend
                                      ip = clist(idom)%lpdx(ipart)
                                      ! check against all particles 
                                      ! in the other cell
                                      DO jpart=jstart,jend
                                          jp = clist(idom)%lpdx(jpart)
                                          dx  = xp(1,ip) - xp(1,jp)
                                          dy  = xp(2,ip) - xp(2,jp)
                                          IF (dimens .GT. 2) THEN
                                              dz  = xp(3,ip) - xp(3,jp)
                                              dij = (dx*dx)+(dy*dy)+(dz*dz)
                                          ELSE
                                              dz = 0.0_MK
                                              dij = (dx*dx)+(dy*dy)
                                          ENDIF
                                          IF (dij .GT. cut2) CYCLE
                                          !-------------------------------------
                                          ! Particle ip interacts with
                                          ! particle jp here... and
                                          ! vice versa to use symmetry.
                                          !-------------------------------------
                                          !
                                          ! USER CODE HERE:
                                          ! dpd(:,ip) = f(pdata(:,ip),
                                          !               pdata(:,jp),dij)
                                          ! dpd(:,jp) = -dpd(:,ip)
                                      ENDDO
                                  ENDDO
                              ENDIF       ! ibox .EQ. jbox
                          ENDDO           ! iinter
                      ENDDO               ! i
                  ENDDO                   ! j
              ENDDO                       ! k
          ENDDO                           ! idom

      !-------------------------------------------------------------------------
      !  PARTICLE-PARTICLE INTERACTIONS not using symmetry
      !-------------------------------------------------------------------------
      ELSE
          DO idom=1,nsublist
              n1  = Nm(1,idom)
              n2  = Nm(1,idom)*Nm(2,idom)
              nz  = Nm(3,idom)
              IF (dimens .EQ. 2) THEN
                  n2 = 0
                  nz = 2
              ENDIF 
              ! loop over all REAL cells (the -2 in the end does this)
              DO k=1,nz-2
                  DO j=1,Nm(2,idom)-2
                      DO i=1,Nm(1,idom)-2
                          ! index of the center box
                          cbox = i + 1 + n1*j + n2*k
                          ! loop over all box-box interactions
                          DO iinter=1,nnp
                              ! determine box indices for this interaction
                              ibox = cbox+(inp(1,iinter)+n1*inp(2,iinter)+ &
     &                               n2*inp(3,iinter))
                              jbox = cbox+(jnp(1,iinter)+n1*jnp(2,iinter)+ &
     &                               n2*jnp(3,iinter))
                              !-------------------------------------------------
                              !  Read indices and check if cell is empty
                              !-------------------------------------------------
                              istart = clist(idom)%lhbx(ibox)
                              iend   = clist(idom)%lhbx(ibox+1)-1
                              IF (iend .LT. istart) CYCLE
                              !-------------------------------------------------
                              !  Do all interactions within the box itself
                              !-------------------------------------------------
                              IF (ibox .EQ. jbox) THEN
                                  DO ipart=istart,iend
                                      ip = clist(idom)%lpdx(ipart)
                                      DO jpart=istart,iend
                                          jp = clist(idom)%lpdx(jpart)
                                          ! No particle interacts with
                                          ! itself
                                          IF (ip .EQ. jp) CYCLE
                                          dx  = xp(1,ip) - xp(1,jp)
                                          dy  = xp(2,ip) - xp(2,jp)
                                          IF (dimens .GT. 2) THEN
                                              dz  = xp(3,ip) - xp(3,jp)
                                              dij = (dx*dx)+(dy*dy)+(dz*dz)
                                          ELSE
                                              dz = 0.0_MK
                                              dij = (dx*dx)+(dy*dy)
                                          ENDIF
                                          IF (dij .GT. cut2) CYCLE
                                          !---------------------------------
                                          ! Particle ip interacts with
                                          ! particle jp.
                                          !---------------------------------
                                          !
                                          ! USER CODE HERE:
                                          ! dpd(:,ip) = f(pdata(:,ip),
                                          !               pdata(:,jp),dij)
                                      ENDDO
                                  ENDDO
                              !-------------------------------------------------
                              !  Do interactions with all neighboring boxes
                              !-------------------------------------------------
                              ELSE
                                  ! get pointers to first and last particle 
                                  jstart = clist(idom)%lhbx(jbox)
                                  jend   = clist(idom)%lhbx(jbox+1)-1
                                  ! skip this iinter if empty
                                  IF (jend .LT. jstart) CYCLE
                                  ! loop over all particles inside this cell
                                  DO ipart=istart,iend
                                      ip = clist(idom)%lpdx(ipart)
                                      ! check against all particles 
                                      ! in the other cell
                                      DO jpart=jstart,jend
                                          jp = clist(idom)%lpdx(jpart)
                                          dx  = xp(1,ip) - xp(1,jp)
                                          dy  = xp(2,ip) - xp(2,jp)
                                          IF (dimens .GT. 2) THEN
                                              dz  = xp(3,ip) - xp(3,jp)
                                              dij = (dx*dx)+(dy*dy)+(dz*dz)
                                          ELSE
                                              dz = 0.0_MK
                                              dij = (dx*dx)+(dy*dy)
                                          ENDIF
                                          IF (dij .GT. cut2) CYCLE
                                          !---------------------------------
                                          ! Particle ip interacts with
                                          ! particle jp.
                                          !---------------------------------
                                          !
                                          ! USER CODE HERE:
                                          ! dpd(:,ip) = f(pdata(:,ip),
                                          !               pdata(:,jp),dij)
                                      ENDDO
                                  ENDDO
                              ENDIF       ! ibox .EQ. jbox
                          ENDDO           ! iinter
                      ENDDO               ! i
                  ENDDO                   ! j
              ENDDO                       ! k
          ENDDO                           ! idom
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_template_comp_pp_cell
