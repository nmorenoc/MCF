#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :           ppm_template_connect_interaction
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine has to be called by the user in order
      !                 to compute the particle connection interactions defined
      !                 by the connections. How to use this template:
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
      !                    (7) set problem dimension, processor rank and
      !                        number of processors in the Initialise
      !                        section.
      !                    (8) change the allocation of the xp2,vp2,... to
      !                        use the correct leading dimension
      !                    (9) change the ifdef USE_MPI to whatever indicates
      !                        the use of MPI in you application
      !                    (10)add bond interactions where needed
      !                        (search for USER CODE)
      !                    (11)Optional: Make argument checking conditional
      !                        on the debug level of your application
      !                    (12)Optional: Make this routine use your own
      !                        error logging, message writing, allocation,
      !                        subroutine start, subroutine stop and
      !                        timing routines.
      !
      !  Input        : id(:)      (I) : local to global particle number mapping
      !                                  This is a must!
      !                 xp(:,:)    (F) : particle co-ordinates
      !                 vp(:,:)    (F) : particle data used for interaction
      !                                  (e.g. velocity)
      !                 !-------------------------------------------------------
      !                 ! USER CODE HERE
      !                 ! mp(:,:)  (F) : particle masses
      !                 ! sp(:,:)  (F) : particle strength
      !                 ! ... and so on, what you need.
      !                 ! Dont forget to modify the argument list below
      !                 ! accordingly
      !                 !-------------------------------------------------------
      !                 Npart      (I) : number of particles on this proc
      !                 cd(:,:)    (I) : connection data
      !                 ldc        (I) : leading dimension of cd(:,:), i.e.
      !                                  the length of the connections
      !                 Ncon       (I) : number of connections on this proc
      !                 lloc       (L) : whether all particles for the local
      !                                  connections are already on this proc:
      !                                    .T. all particles are local
      !                                    .F. some particles are on other procs
      !
      !  Input/output : fp(:,:)    (F) : Change of particle data due
      !                                  to interaction.
      !                 !-------------------------------------------------------
      !                 ! USER CODE HERE
      !                 ! Add any output arrays you need.
      !                 ! Dont forget to modify the argument list below
      !                 ! accordingly
      !                 !-------------------------------------------------------
      !                 params(:)  (F) : user defined parameters and/or
      !                                  output from the interaction (i.e.
      !                                  potential energy)
      !  Output       : info       (I) : return status
      !
      !  Routines     : ppm_map_part_ring_shift
      !                 ppm_map_part
      !
      !  Remarks      : Search for USER in this file and add as
      !                 many arguments used for particle interaction to the
      !                 argument list you need for your specific computation.
      !
      !                 The actual form is an example for an interaction that
      !                 takes the particle positions and velocities and gives
      !                 back the resulting forces.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Revision: 1.7 $
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Add as many arguments to the list as you need.
      !
      !  SUBROUTINE ppm_template_connect_interaction(id,xp,vp,Npart,cd,lda,   &
      !                                              Ncon,lloc,fp,params,info)
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_template_connect_interaction(id,xp,vp,Npart,cd,ldc,Ncon, &
     &                                            lloc,fp,params,info)

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      ! USER: USE statements here...
      USE ppm_module_map
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
      INTEGER , DIMENSION(:)  , INTENT(IN   ) :: id
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: vp
      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Add your input arguments here, e.g.
      !  REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: mp
      !  REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: sp
      !  ...
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: Npart
      INTEGER , DIMENSION(:,:), INTENT(IN   ) :: cd
      INTEGER                 , INTENT(IN   ) :: ldc
      INTEGER                 , INTENT(IN   ) :: Ncon
      LOGICAL                 , INTENT(IN   ) :: lloc
      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Add your output arguments here
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(INOUT) :: fp
      REAL(MK), DIMENSION(:)  , INTENT(INOUT) :: params
      INTEGER                 , INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                           :: i,j,k,l,m,dims,pcount,psize,found
      INTEGER                           :: itarget,isource,hops,Lpart,nproc,rank
      INTEGER                           :: lda
      INTEGER , DIMENSION(ldc)          :: p
      INTEGER , DIMENSION(:)  , POINTER :: id2
      INTEGER , DIMENSION(:,:), POINTER :: cd_local,pmap
      REAL(MK), DIMENSION(:,:), POINTER :: xp_con,vp_con,fp_con
      ! USER: if allocation of the following varaibles fails due to stack 
      ! size limitations, try putting them in a module and add a USE
      ! statement for it above.
      REAL(MK), DIMENSION(:,:), POINTER :: xp2,vp2
      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Add the variable names for the copies here, e.g.
      !  REAL(MK), DIMENSION(:,:), POINTER :: mp2,sp2
      !  ...
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !  USER: Set problem dimension, processor rank and number of
      !  processors here...
      !-------------------------------------------------------------------------
      info = 0
      dims = SIZE(xp,1)
      lda  = SIZE(fp,1)
      rank =
      nproc =

      !-------------------------------------------------------------------------
      !  Check input arguments.
      !  USER: you may want to make this conditional on the debug level of
      !  your application...
      !-------------------------------------------------------------------------
      IF (Ncon .LT. 0) THEN
          WRITE(*,'(2A)') '(ppm_template_connect_interaction): ', &
     &        'Ncon must be >=0'
          info = -1
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Allocate memory for the copy of the local particles
      !-------------------------------------------------------------------------
      ALLOCATE(id2(Npart),xp2(dims,Npart),vp2(dims,Npart),STAT=info)
      IF (info .NE. 0) THEN
          WRITE(*,'(2A,I8)') '(ppm_template_connect_interaction): ',  &
     &        'allocation failed on line ',__LINE__
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Make a copy of all local particles including all additional infos
      !  that they carry
      !-------------------------------------------------------------------------
      Lpart = Npart
      id2(:)   = id(:)
      xp2(:,:) = xp(:,:)
      vp2(:,:) = vp(:,:)
      !  mp2(:,:) = mp(:,:)
      !  sp2(:,:) = sp(:,:)
      !  ...

      !-------------------------------------------------------------------------
      !  Count how many different particles are involved in the connections
      !-------------------------------------------------------------------------
      ALLOCATE(cd_local(ldc,Ncon),STAT=info)
      IF (info .NE. 0) THEN
          WRITE(*,'(2A,I8)') '(ppm_template_connect_interaction): ',  &
     &        'allocation failed on line ',__LINE__
          GOTO 9999
      ENDIF

      cd_local(:,:) = cd(:,:)

      pcount = 0
      DO j = 1,Ncon
         DO i = 1,ldc
            IF (cd_local(i,j) .NE. -1) THEN
                pcount = pcount + 1
                found = cd_local(i,j)
                WHERE(cd_local .EQ. found) cd_local = -1
            ENDIF
         ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  USER CODE HERE:
      !  Allocate memory for the particles that are involved in the connections
      !-------------------------------------------------------------------------
      ALLOCATE(pmap(3,pcount),xp_con(dims,pcount),vp_con(dims,pcount), &
     &         fp_con(lda,pcount),STAT=info)
      IF (info .NE. 0) THEN
          WRITE(*,'(2A,I8)') '(ppm_template_connect_interaction): ',       &
     &        'allocation failed on line ',__LINE__
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Gather all particles that are inside a local connection (local)
      !-------------------------------------------------------------------------
      pcount = 0
      DO j = 1,Ncon
         DO i = 1,ldc
            DO l = 1,Npart
               IF (cd(i,j) .EQ. id(l)) THEN
                   !------------------------------------------------------------
                   !  Check if we already have this particle
                   !------------------------------------------------------------
                   found = 0
                   DO m = 1,pcount
                      IF (pmap(3,m) .EQ. cd(i,j)) THEN
                          found = 1
                          EXIT
                      ENDIF
                   ENDDO

                   !------------------------------------------------------------
                   !  If yes, go to the next one
                   !------------------------------------------------------------
                   IF (found .EQ. 1) CYCLE

                   pcount = pcount + 1

                   !------------------------------------------------------------
                   !  DO NOT EDIT
                   !  Store mapping information in order to be able to send the
                   !  results back where they belong to.
                   !------------------------------------------------------------
                   pmap(1,pcount) = rank
                   pmap(2,pcount) = l
                   pmap(3,pcount) = cd(i,j)

                   !------------------------------------------------------------
                   !  USER CODE HERE:
                   !  Copy the particle data
                   !------------------------------------------------------------
                   xp_con(:,pcount) = xp(:,l)
                   vp_con(:,pcount) = vp(:,l)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

#ifdef USE_MPI
      !-------------------------------------------------------------------------
      !  Send the particles around the ring in order to get the missing ones
      !  for the local connections
      !-------------------------------------------------------------------------
      IF (.NOT.lloc) THEN
          !---------------------------------------------------------------------
          !  Compute whom to send the connections and who to receive the
          !  connections from
          !---------------------------------------------------------------------
          itarget = MOD(nproc + rank - 1, nproc)
          isource = MOD(rank + 1, nproc)

          !---------------------------------------------------------------------
          !  Send it to all processors
          !---------------------------------------------------------------------
          hops = nproc - 1

          !---------------------------------------------------------------------
          !  Send the particle data through the ring
          !---------------------------------------------------------------------
          DO k = 1,hops
             !------------------------------------------------------------------
             !  USER CODE HERE
             !  Shift the particle data to the neighbour processor and receive
             !  it from the other neighbour. Add any additional particle data
             !  you need. id2 needs to be send around too!
             !------------------------------------------------------------------
             CALL ppm_map_part_ring_shift(xp2,dims,Lpart,itarget,isource,info)
             CALL ppm_map_part(vp2,dims,Lpart,Lpart,-1,ppm_param_map_push,info)
             CALL ppm_map_part(id2,Lpart,Lpart,-1,ppm_param_map_push,info)
             CALL ppm_map_part(xp2,dims,Lpart,Lpart,-1,ppm_param_map_send,info)
             CALL ppm_map_part(id2,Lpart,Lpart,-1,ppm_param_map_pop,info)
             CALL ppm_map_part(vp2,dims,Lpart,Lpart,-1,ppm_param_map_pop,info)
             CALL ppm_map_part(xp2,dims,Lpart,Lpart,-1,ppm_param_map_pop,info)

             !------------------------------------------------------------------
             ! Gather all particles that are inside a local connection (copies)
             !------------------------------------------------------------------
             DO j = 1,Ncon
                DO i = 1,ldc
                   DO l = 1,Lpart
                      IF (cd(i,j) .EQ. id2(l)) THEN
                          !-----------------------------------------------------
                          !  Check if we already have this particle
                          !-----------------------------------------------------
                          found = 0
                          DO m = 1,pcount
                             IF (pmap(3,m) .EQ. cd(i,j)) THEN
                                 found = 1
                                 EXIT
                             ENDIF
                          ENDDO

                          !-----------------------------------------------------
                          !  If yes, go to the next one
                          !-----------------------------------------------------
                          IF (found .EQ. 1) CYCLE

                          pcount = pcount + 1

                          !-----------------------------------------------------
                          !  DO NOT EDIT
                          !  Store mapping information in order to be able to
                          !  send the results back where they belong to.
                          !-----------------------------------------------------
                          pmap(1,pcount) = MOD(rank + k,nproc)
                          pmap(2,pcount) = l
                          pmap(3,pcount) = cd(i,j)

                         !------------------------------------------------------
                         !  USER CODE HERE:
                         !  Copy the particle data
                         !------------------------------------------------------
                         xp_con(:,pcount) = xp2(:,l)
                         vp_con(:,pcount) = vp2(:,l)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
         ENDDO
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  Translate the connection data (cd) from a global particle numbering
      !  to a local numbering (cd_local). This makes life easier for computing
      !  the interactions
      !-------------------------------------------------------------------------
      DO j = 1,Ncon
         DO i = 1,ldc
            DO k = 1,pcount
               IF (pmap(3,k) .EQ. cd(i,j)) THEN
                   cd_local(i,j) = k
                   EXIT
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      !------------------------------------------------------------------------
      !  USER CODE HERE
      !  Loop over all connections
      !
      !  cd_local(1:ldc,1:Ncon)
      !  ldc  : the length of the connection (i.e. number of particles in this
      !         connection)
      !  Ncon : the number of connections
      !
      !  e.g. lets take connection number 2
      !    p(:) = cd_local(:,2)
      !  and a connection length of 3, then p contains the local number of
      !  the particles in this connection:
      !  xp_con(:,p(1)) is the first particle in the connection
      !  xp_con(:,p(2)) is the second particle in the connection
      !  xp_con(:,p(3)) is the third particle in the connection
      !
      !------------------------------------------------------------------------
      fp_con(:,:) = 0.0_MK

      DO j = 1,Ncon
         p(:) = cd_local(:,j)

         !DO i = 2,ldc
         !    dx = xp_con(1,p(i-1)) - xp_con(1,p(i))
         !    ...
         !    fp_con(1,p(i-1)) = fp_con(1,p(i-1)) + dx*fc
         !    ...
         !ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Add up the results to the existing ones (only the ones we have the
      !  particles for)
      !-------------------------------------------------------------------------
      DO j = 1,pcount
         IF (pmap(1,j) .EQ. rank) THEN
             DO i = 1,lda
                fp(i,pmap(2,j)) = fp(i,pmap(2,j)) + fp_con(i,j)
             ENDDO
         ENDIF
      ENDDO

#ifdef USE_MPI
      !-------------------------------------------------------------------------
      !  Send the results back where they belong to only if requested
      !-------------------------------------------------------------------------
      IF (.NOT.lloc) THEN
          !---------------------------------------------------------------------
          !  USER CODE HERE
          !  Send the forces or whatever back where it belongs to and add it to
          !  the existing ones. pmap _must_ to be send along too because it
          !  contains the mapping informations
          !---------------------------------------------------------------------
          DO k = 1,hops
             CALL ppm_map_part_ring_shift(fp_con,lda,pcount,itarget,isource,info)
             CALL ppm_map_part(pmap,3,pcount,pcount,-1,ppm_param_map_push,info)
             CALL ppm_map_part(fp_con,lda,pcount,pcount,-1,ppm_param_map_send,info)
             CALL ppm_map_part(pmap,3,pcount,pcount,-1,ppm_param_map_pop,info)
             CALL ppm_map_part(fp_con,lda,pcount,pcount,-1,ppm_param_map_pop,info)

             DO j = 1,pcount
                IF (pmap(1,j) .EQ. rank) THEN
                    DO i = 1,lda
                       fp(i,pmap(2,j)) = fp(i,pmap(2,j)) + fp_con(i,j)
                    ENDDO
                ENDIF
             ENDDO
          ENDDO
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Deallocate memory or you will have a serious leak
      !-------------------------------------------------------------------------
      DEALLOCATE(cd_local,pmap,fp_con,vp_con,xp_con,xp2,vp2,id2,STAT=info)
      IF (info .NE. 0) THEN
          WRITE(*,'(2A,I8)') '(ppm_template_connect_interaction): ',  &
     &        'deallocation failed on line ',__LINE__
          GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE

      RETURN

      END SUBROUTINE ppm_template_connect_interaction

