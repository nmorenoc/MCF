#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :              ppm_template_ring_interaction
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine has to be called by the user in order
      !                 to compute the N**2 particle interactions on a
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
      !                    (7) set problem dimensions, processor rank and
      !                        number of processors in the Initialise
      !                        section.
      !                    (8) change the allocation of the xp2,vp2,... to
      !                        use the correct leading dimension
      !                    (9) change the ifdef USE_MPI to whatever indicates
      !                        the use of MPI in you application
      !                    (10)add particle-particle interactions where
      !                        needed (search for USER CODE)
      !                    (11)replace all CALLs to
      !                        ppm_template_comp_pp_ring with your PP
      !                        interaction subroutine
      !                    (12)Optional: Make argument checking conditional
      !                        on the debug level of your application
      !                    (13)Optional: Make this routine use your own
      !                        error logging, message writing, allocation,
      !                        subroutine start, subroutine stop and
      !                        timing routines.
      !
      !  Input        : xp(:,:)     (F) : particle co-ordinates
      !                 vp(:,:)     (F) : particle data used for interaction
      !                                   (e.g. velocity)
      !                 !-------------------------------------------------------
      !                 ! USER CODE HERE
      !                 ! mp(:,:)   (F) : particle masses
      !                 ! sp(:,:)   (F) : particle strength
      !                 ! ... and so on, what you need.
      !                 ! Dont forget to modify the argument list below
      !                 ! accordingly
      !                 !-------------------------------------------------------
      !                 Npart       (I) : number of particles on this proc
      !                 lsymm       (L) : Whether to use symmetry or not:
      !                                    .F.  pp interaction w/o symmetry
      !                                    .T.  pp interaction w/  symmetry
      !
      !  Input/output : params(:)   (F) : user defined parameters and/or
      !                                   output from the interaction (i.e.
      !                                   potential energy)
      !
      !  Output       : fp(:,:)     (F) : Change of particle data (vp) due
      !                                   to interaction.
      !                 !-------------------------------------------------------
      !                 ! USER CODE HERE
      !                 ! Add any output arrays you need.
      !                 ! Dont forget to modify the argument list below
      !                 ! accordingly
      !                 !-------------------------------------------------------
      !                 info        (I) : return status
      !
      !  Routines     : ppm_template_comp_pp_ring
      !                 ppm_map_part_ring_shift
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
      !  $Revision: 1.17 $
      !  Revision 1.3  2004/04/23 17:23:45  oingo
      !  Revision 1.2  2004/04/22 11:40:36  oingo
      !  Revision 1.1  2004/04/22 08:24:43  oingo
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Add as many arguments to the list as you need. Be aware that the second
      !  dimension of the arrays (i.e. UBOUND(array,2)) has to be Npart, e.g.
      !
      !  SUBROUTINE ppm_template_ring_interaction(xp,vp,mp,sp,lsymm,fp,params,&
      ! &                                           info)
      !-------------------------------------------------------------------------
      SUBROUTINE ppm_template_ring_interaction(xp,vp,Npart,lsymm,fp,params, &
     &                                         info)

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
      LOGICAL                 , INTENT(IN   ) :: lsymm
      REAL(MK), DIMENSION(:,:), INTENT(  OUT) :: fp
      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Add your output arguments here
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:)  , INTENT(INOUT) :: params
      INTEGER                 , INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER                           :: i,j,dims,lda,rank,nproc
      INTEGER                           :: hops,isource,itarget
      INTEGER                           :: ll,lu,rl,ru
      INTEGER                           :: Lpart
      ! USER: if allocation of the following three varaibles fails due to stack 
      ! size limitations, try putting them in a module and add a USE
      ! statement for it above.
      REAL(MK), DIMENSION(:,:), POINTER :: xp2,vp2,fp2
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
      fp = 0.0_MK
      dims = SIZE(xp,1)           ! space dimension
      lda  = SIZE(fp,1)           ! leading dimension of particle data
      rank =
      nproc =

      !-------------------------------------------------------------------------
      !  Check input arguments.
      !  USER: you may want to make this conditional on the debug level of
      !  your application...
      !-------------------------------------------------------------------------
      IF (Npart .LE. 0) THEN
         WRITE(*,'(2A)') '(ppm_template_ring_interaction): ', &
     &        'Npart must be >0'
         info = -1
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Allocate memory for the copy of the local particles
      !-------------------------------------------------------------------------
      ALLOCATE(xp2(dims,Npart),vp2(lda,Npart),fp2(lda,Npart),STAT=info)
      IF (info .NE. 0) THEN
         WRITE(*,'(2A,I8)') '(ppm_template_ring_interaction): ',  &
     &        'allocation failed on line ',__LINE__
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Make a copy of all local particles including all additional infos
      !  that they carry
      !-------------------------------------------------------------------------
      Lpart = Npart
      xp2(:,:) = xp(:,:)
      vp2(:,:) = vp(:,:)
      fp2(:,:) = fp(:,:)
      !  mp2(:,:) = mp(:,:)
      !  sp2(:,:) = sp(:,:)
      !  ...

      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Calculate the interactions of the local particles with themselves.
      !
      !  Here you have to adjust the name of the subroutine (if you changed it)
      !  and the arguments in the list according to the data you need.
      !-------------------------------------------------------------------------
      CALL ppm_template_comp_pp_ring(xp ,vp ,fp ,Npart,          &
     &                               xp2,vp2,fp2,Lpart,          &
     &                               lsymm,params,1,info)
      IF (info .NE. 0) THEN
         WRITE(*,'(2A,I8)') '(ppm_template_ring_interaction):',   &
     &        'interaction calculation failed on line ',__LINE__
         GOTO 9999
      ENDIF

#ifdef USE_MPI
      !-------------------------------------------------------------------------
      !  Compute whom to send the copy and who to receive the copy from
      !-------------------------------------------------------------------------
      itarget = MOD(nproc + rank - 1, nproc)
      isource = MOD(rank + 1, nproc)
      
      !-------------------------------------------------------------------------
      !  Compute how often we have to shift a copy of the particles to the
      !  neighbour processor
      !-------------------------------------------------------------------------
      IF (lsymm) THEN
         IF (MOD(nproc,2) .EQ. 0) THEN
            hops = (nproc / 2) - 1
         ELSE
            hops = ((nproc + 1) / 2) - 1
         ENDIF
      ELSE
         hops = nproc - 1
      ENDIF

      DO i = 1,hops
         CALL ppm_map_part_ring_shift(xp2,dims,Lpart,itarget,isource,info)
         !----------------------------------------------------------------------
         !  USER CODE HERE
         !  Push any array you want to be sent along with the particle on the
         !  stack, but only the copy-arrays. The particle positions dont need
         !  to be pushed once more since theyve already been pushed on the
         !  stack by the ppm_map_part_ring_shift subroutine.
         !  After pushing everything, send and retrieve the arrays of the other
         !  processor by calling the pop subroutine in reverse order. xp2 has
         !  to be poped explicitly.
         !  Dont forget to modify the call of the interaction calculation
         !  subroutine as above.
         !----------------------------------------------------------------------
         CALL ppm_map_part(vp2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
         CALL ppm_map_part(fp2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
         CALL ppm_map_part(xp2,dims,Lpart,Lpart,-1,ppm_param_map_send,info)
         CALL ppm_map_part(fp2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
         CALL ppm_map_part(vp2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
         CALL ppm_map_part(xp2,dims,Lpart,Lpart,-1,ppm_param_map_pop,info)

         CALL ppm_template_comp_pp_ring(xp ,vp ,fp ,Npart,          &
     &                                  xp2,vp2,fp2,Lpart,          &
     &                                  lsymm,params,0,info)
         IF (info .NE. 0) THEN
            WRITE(*,'(2A,I8)') '(ppm_template_ring_interaction): ',  &
     &           'interaction calculation failed on line ',__LINE__
            GOTO 9999
         ENDIF
      ENDDO
      
      !-------------------------------------------------------------------------
      !  If we have symmetry we only have to send the particles half the round,
      !  so we have to check the case of an even number of processors
      !-------------------------------------------------------------------------
      IF (lsymm .AND. (MOD(nproc,2) .EQ. 0)) THEN
         CALL ppm_map_part_ring_shift(xp2,dims,Lpart,itarget,isource,info)
         !---------------------------------------------------------------------
         !  USER CODE HERE
         !  Do the push-send-pop business as above
         !---------------------------------------------------------------------
         CALL ppm_map_part(vp2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
         CALL ppm_map_part(fp2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
         CALL ppm_map_part(xp2,dims,Lpart,Lpart,-1,ppm_param_map_send,info)
         CALL ppm_map_part(fp2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
         CALL ppm_map_part(vp2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
         CALL ppm_map_part(xp2,dims,Lpart,Lpart,-1,ppm_param_map_pop,info)

         !---------------------------------------------------------------------
         !  The processor with the higher rank computes the upper half of
         !  the particles and the opposite processor the lower half
         !---------------------------------------------------------------------
         IF (rank .GT. hops) THEN
            ll = (Npart / 2) + 1
            lu = Npart
            rl = 1
            ru = Lpart
         ELSE
            ll = 1
            lu = Npart
            rl = 1
            ru = Lpart / 2
         ENDIF

         CALL ppm_template_comp_pp_ring(xp(:,ll:lu),vp(:,ll:lu),fp(:,ll:lu),  &
     &                                  lu-ll+1,xp2(:,rl:ru),vp2(:,rl:ru),    &
     &                                  fp2(:,rl:ru),ru-rl+1,lsymm,params,    &
     &                                  0,info)
         IF (info .NE. 0) THEN
            WRITE(*,'(2A,I8)') '(ppm_template_ring_interaction): ',&
     &           'interaction calculation failed on line ',__LINE__
            GOTO 9999
         ENDIF
      ENDIF
      
      !-------------------------------------------------------------------------
      !  Send the particles where they belong to
      !-------------------------------------------------------------------------
      IF (lsymm) THEN
         itarget = MOD(rank + (nproc / 2),nproc)
         isource = MOD(nproc + rank - (nproc / 2),nproc)
      ENDIF

      IF (itarget .NE. rank) THEN
         CALL ppm_map_part_ring_shift(xp2,dims,Lpart,itarget,isource,info)
         !---------------------------------------------------------------------
         !  USER CODE HERE
         !  Do the push-send-pop business as above
         !---------------------------------------------------------------------
         CALL ppm_map_part(vp2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
         CALL ppm_map_part(fp2,lda,Lpart,Lpart,-1,ppm_param_map_push,info)
         CALL ppm_map_part(xp2,dims,Lpart,Lpart,-1,ppm_param_map_send,info)
         CALL ppm_map_part(fp2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
         CALL ppm_map_part(vp2,lda,Lpart,Lpart,-1,ppm_param_map_pop,info)
         CALL ppm_map_part(xp2,dims,Lpart,Lpart,-1,ppm_param_map_pop,info)

         IF (Lpart .NE. Npart) THEN
            info = -2
            WRITE(*,'(2A,I8)') '(ppm_template_ring_interaction): ',  &
     &           'Not all particles came back on line ',__LINE__
            GOTO 9999
         ENDIF
      ENDIF
#endif

      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Add the particle changes of the local and the long traveled copy
      !-------------------------------------------------------------------------
      DO j = 1,lda
         DO i = 1,Npart
            fp(j,i) = fp(j,i) + fp2(j,i)
         ENDDO
      ENDDO
      
      !-------------------------------------------------------------------------
      !  USER CODE HERE
      !  Deallocate the memory of the copies or you will have a serious leak
      !-------------------------------------------------------------------------
      DEALLOCATE(fp2,vp2,xp2,STAT=info)
      IF (info .NE. 0) THEN
         WRITE(*,'(2A,I8)') '(ppm_template_ring_interaction): ',  &
     &        'deallocation failed on line ',__LINE__
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      RETURN
      END SUBROUTINE ppm_template_ring_interaction
