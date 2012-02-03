#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_solve_3d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Poisson solver  using FFTW
      !                 Solves the negative poisson equation in a 3-dimensional 
      !                 periodic domain:
      !                          - Laplacian of Phi = omega
      !                 It takes the field quantity omega in DATA_fv which is 
      !                 assumed to be periodic outside the domain. The most 
      !                 efficient initial topology for the fieldsolver is 
      !                 a x-pencil topology. After performing a FFT on the  
      !                 x-pencils the data is mapped onto y-pencils and 
      !                 later onto z-pencils.
      !                 The version solve_init is based on recomputed FFT-plans
      !                 created in ppm_fdsolver_init and destroyed in
      !                 ppm_fdsolver_finalize.
      !                 In the slab version (solve_slab) the FFTs are performed
      !                 on xy slabs and on z-penzils to save field mappings. 
      !                 This version requires the call of ppm_fdsolver_init and
      !                 ppm_fdsolver_finalize.
      !                 The poisson equation is solved in the Fourier space 
      !                 and the result transformed backward.
      !                 The solution Phi is finally returned in DATA_fv.
      !                 Note: field quantity must live on current topology
      !
      !                 Usage:   
      !    
      !                    ppm_fdsolver_init(arguments)
      !   
      !                    ppm_fdsolver_solve_init(arguments)
      !                 or ppm_fdsolver_solve_slab(arguments)
      !   
      !                    ppm_fdsolver_finalize(info)
      !   
      !   
      !                 Standalone Usage:
      !
      !                    ppm_fdsolver_solve(arguments)
      !
      !
      !
      !
      !  Input        : 
      !                  lda_fv      (I) size of leading dimension in VectorCase
      !                  mesh_id_user(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) first topology in
      !                                               solve:      xpencils  
      !                                               solve_init: xpencils  
      !                                               solve_slab: xy slab  
      !                                  topo_ids(2) second  topology in
      !                                               solve:      ypencils  
      !                                               solve_init: ypencils  
      !                                               solve_slab: zpencils  
      !                                  topo_ids(3) third   topology in
      !                                               solve:      zpencils  
      !                                               solve_init: xpencils  
      !                                               solve_slab: not used  
      !
      !                  mesh_ids(3) (I) mesh IDs where the FFTs are performed
      !                                  mesh_ids(1) first mesh 
      !                                             solve:      xpencils,real  
      !                                             solve_init: xpencils,real  
      !                                             solve_slab: xy slab ,real 
      !                                  mesh_ids(2) second mesh 
      !                                             solve:      xpencils,complex
      !                                             solve_init: xpencils,complex  
      !                                             solve_slab: xy slab ,complex
      !                                  mesh_ids(3) third mesh 
      !                                             solve:      ypencils,complex
      !                                             solve_init: ypencils,complex  
      !                                             solve_slab: zpencils,complex
      !                                  mesh_ids(4) forth mesh 
      !                                             solve:      zpencils,complex
      !                                             solve_init: zpencils,complex  
      !                                             solve_slab: not used
      !                  ghostsize(3) (I)ghostsize
      !                                
      !
      !  Input/output : 
      !                  DATA_fv(:,:,:,:) (F) field data         
      !                                
      !                                
      !
      !  Output       : info       (I) return status. =0 if no error.
      !                   
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_solve_3d.f,v $
      !  Revision 1.28  2006/09/04 18:34:45  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.26  2006/04/10 08:54:30  pchatela
      !  Made xp a REAL, DIMENSION(:,:), POINTER to get rid of warnings
      !
      !  Revision 1.25  2006/04/07 17:41:29  hiebers
      !  Changed type of variable xp to POINTER
      !
      !  Revision 1.24  2005/08/23 15:23:12  hiebers
      !  Bugfix in scalar version by adjusting an argument in the call 
      ! of the FFT       routines
      !
      !  Revision 1.23  2005/08/03 14:35:12  ivos
      !  Shortened the subrountine names to meet the F90 31-character limit
      !  standard. pgf90 on gonzales had problems...
      !
      !  Revision 1.22  2005/02/16 12:38:20  hiebers
      !  bugfix in preprocessor if statement
      !
      !  Revision 1.21  2005/02/16 12:01:07  hiebers
      !  Major addition: finalized scalar and vector version,
      !        implemented version <solve_init> that uses FFT plans created
      !        in ppm_fdsolver_init and destroyed in ppm_fdsolver_finalize
      !        whereas the version <solve> remains a standalone solver,
      !        implemented version <solve_slab> that uses xy-slab topology
      !        and z-pencil topology (needs ppm_fdsolver_init and finalize)
      !
      !  Revision 1.20  2005/02/01 16:38:56  hiebers
      !  reduced size of DATA_fv_com
      !
      !  Revision 1.19  2004/11/03 11:14:55  hiebers
      !  Exchanged __SXF90 by __MATHKEISAN
      !
      !  Revision 1.18  2004/10/01 16:08:59  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.17  2004/08/27 09:45:15  hiebers
      !  added ghostsize as an argument
      !
      !  Revision 1.16  2004/08/26 10:45:52  hiebers
      !  eliminated correction of the field margin. correction is now done in
      !  the fft routines.
      !
      !  Revision 1.9  2004/08/19 13:14:28  hiebers
      !  debugged scalar/vector version
      !
      !  Revision 1.8  2004/08/18 16:32:43  hiebers
      !  intermediate status of debugging vector version
      !
      !  Revision 1.7  2004/08/18 13:35:33  hiebers
      !  bug fix for non xpencil-topology version
      !
      !  Revision 1.3  2004/07/26 15:38:47  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.2  2004/07/26 11:59:38  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.1  2004/07/26 08:52:29  hiebers
      !  Recommited, formerly ppm_module_fieldsolver
      !
      !  Revision 1.1  2004/05/19 15:32:46  hiebers
      !  implementation from scratch
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __CASE == __SLAB

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_slab_3d_ss(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_slab_3d_sd(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize,info)
#endif

#endif


#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_slab_3d_vs(DATA_fv,lda_fv,mesh_id_user, &
        topo_ids,mesh_ids,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_slab_3d_vd(DATA_fv,lda_fv,mesh_id_user, & 
        topo_ids,mesh_ids,ghostsize,info)
#endif
#endif


#elif   __CASE == __INIT

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_init_3d_ss(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_init_3d_sd(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize,info)
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_init_3d_vs(DATA_fv,lda_fv,mesh_id_user, &
        topo_ids,mesh_ids,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_init_3d_vd(DATA_fv,lda_fv,mesh_id_user, & 
        topo_ids,mesh_ids,ghostsize,info)
#endif
#endif



#else
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_3d_ss(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_3d_sd(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize,info)
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_3d_vs(DATA_fv,lda_fv,mesh_id_user, &
        topo_ids,mesh_ids,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_3d_vd(DATA_fv,lda_fv,mesh_id_user, & 
        topo_ids,mesh_ids,ghostsize,info)
#endif
#endif

#endif


      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"


      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------

      USE   ppm_module_data
      USE   ppm_module_data_mesh
      USE   ppm_module_substart
      USE   ppm_module_substop
      USE   ppm_module_write
      USE   ppm_module_error
      USE   ppm_module_alloc
      USE   ppm_module_mktopo
      USE   ppm_module_mesh_define
      USE   ppm_module_fdsolver_map
      USE   ppm_module_util_fft_forward
      USE   ppm_module_util_fft_backward
      USE   ppm_module_fdsolver_poisson
      USE   ppm_module_fdsolver_fft_fd
      USE   ppm_module_fdsolver_fft_bd
 

      IMPLICIT NONE


#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! data
#if   __DIM == __SFIELD
      REAL(MK), DIMENSION(:,:,:,:),    POINTER     :: DATA_fv
#elif __DIM == __VFIELD
      REAL(MK), DIMENSION(:,:,:,:,:),  POINTER     :: DATA_fv
#endif
      ! mesh ID of the data
      INTEGER                    , INTENT(IN)      :: mesh_id_user

#if   __DIM == __VFIELD
      INTEGER                    , INTENT(IN)      :: lda_fv
#endif

      ! topo / mesh ID for the mapping
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(4)      , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: ghostsize
      INTEGER                    , INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                   :: t0
      ! parameter for alloc
      INTEGER, DIMENSION(3)                      :: lda
      ! counters
      INTEGER                                    :: i,j,k,l,m 
      INTEGER                                    :: ires, jres, lres
      ! 1/number of gridpoints
      REAL(MK)                                   :: rN
      ! result
 
#if   __DIM == __SFIELD
      INTEGER, PARAMETER                         :: lda_fv = 1
#endif

#if   __DIM == __SFIELD
      INTEGER , DIMENSION(4  )                   :: lda_DATA_fv
      COMPLEX(MK), DIMENSION(:,:,:,:), POINTER   :: DATA_fv_com
      INTEGER , DIMENSION(4  )                   :: lda_DATA_fv_com
#elif __DIM == __VFIELD
      INTEGER , DIMENSION(5  )                   :: lda_DATA_fv
      COMPLEX(MK), DIMENSION(:,:,:,:,:), POINTER :: DATA_fv_com
      INTEGER , DIMENSION(5  )                   :: lda_DATA_fv_com
#endif

      INTEGER, DIMENSION(2)                      :: topo_ids_tmp
      INTEGER, DIMENSION(2)                      :: mesh_ids_tmp

      REAL(MK), DIMENSION(:,:,:), POINTER        :: data_in
      COMPLEX(MK), DIMENSION(:,:,:), POINTER     :: data_com
      COMPLEX(MK), DIMENSION(:,:,:), POINTER     :: FFT_x, FFT_xy, FFT_xyz
      REAL(MK), DIMENSION(:,:,:), POINTER        :: Result


      ! variables
      REAL(MK), DIMENSION(:,:),POINTER        :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp, asign
      REAL(MK), DIMENSION(3  )                :: min_phys, max_phys
      REAL(MK), DIMENSION(3  )                :: length
      REAL(MK), DIMENSION(3  )                :: length_phys
      INTEGER , DIMENSION(6  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id, mesh_id
      INTEGER                                 :: yhmax, zhmax
      INTEGER                                 :: mesh_id_internal
      INTEGER                                 :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      INTEGER                                 :: mesh_id_zpen
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart, istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen, istart_trans
      INTEGER , DIMENSION(:,:), POINTER       :: istart_zpen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata, ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen, ndata_trans
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_zpen
      INTEGER , DIMENSION(:  ), POINTER       :: sub2proc
      INTEGER , DIMENSION(:  ), POINTER       :: isublist
      INTEGER                                 :: nsublist, idom
      INTEGER                                 :: dim, n
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_user, topo_id_internal
      INTEGER                                 :: topo_id_xpen, topo_id_ypen
      INTEGER                                 :: topo_id_zpen
      INTEGER, DIMENSION(3)                   :: Nm, Nm_com, Nm_poisson
      LOGICAL                                 :: Its_xpencil_topo
      LOGICAL                                 :: Its_xyslab_topo
      CHARACTER(LEN=ppm_char)                 :: mesg

      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_solve_3d',t0,info)

      !-------------------------------------------------------------------------
      ! Check arguments
      !-------------------------------------------------------------------------
      
      IF (ppm_debug .GT. 0) THEN
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_fdsolver_solve',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ppm_field_topoid .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_solve',  &
     &            'No field topology has been defined so far',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF


#if  !(defined(HAVE_LIBFFTW3) | defined(__MATHKEISAN))

      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error


#ifndef HAVE_LIBFFTW3
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_solve_3d',  &
     &    'PPM was compiled without fftw support',__LINE__,info)
#endif

#ifndef __MATHKEISAN
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_solve_3d',  &
     &    'PPM was compiled without MATHKEISAN support',__LINE__,info)
#endif

      GOTO 9999   


#else


      !----------------------------------------------------------------------
      !  Allocate the isublist array
      !----------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      lda(1) = ppm_nsublist(ppm_field_topoid)
      CALL ppm_alloc(isublist,lda,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         WRITE (mesg,'(A,I10,A)') 'allocating ',      &
     &                   ppm_nsublist(ppm_field_topoid) ,' isublist failed'
         CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve_3d',      &
     &           mesg,__LINE__,info)
         GOTO 9999
      ENDIF
         


      !-------------------------------------------------------------------------
      !  Initialize variables
      !-------------------------------------------------------------------------

      bcdef(1)          = ppm_param_bcdef_periodic
      bcdef(2)          = ppm_param_bcdef_periodic
      bcdef(3)          = ppm_param_bcdef_periodic 
      bcdef(4)          = ppm_param_bcdef_periodic
      bcdef(5)          = ppm_param_bcdef_periodic
      bcdef(6)          = ppm_param_bcdef_periodic

      asign            = ppm_param_assign_internal
      topo_id_internal  = ppm_field_topoid
      topo_id_user      = ppm_user_topoid(topo_id_internal)
      mesh_id_internal  = ppm_meshid(topo_id_internal)%internal(mesh_id_user)
      nsublist          = ppm_nsublist(topo_id_internal)
      isublist          = ppm_isublist(:,topo_id_internal)
      Nm(1)             = ppm_cart_mesh(mesh_id_internal,topo_id_internal)%Nm(1)
      Nm(2)             = ppm_cart_mesh(mesh_id_internal,topo_id_internal)%Nm(2)
      Nm(3)             = ppm_cart_mesh(mesh_id_internal,topo_id_internal)%Nm(3)
      Nm_com(1)         = Nm(1)/2+1
      Nm_com(2)         = Nm(2)
      Nm_com(3)         = Nm(3)
      Npart             = 0

#if   __KIND == __SINGLE_PRECISION
      
      min_phys(1) = ppm_min_physs(1, topo_id_internal) 
      min_phys(2) = ppm_min_physs(2, topo_id_internal) 
      min_phys(3) = ppm_min_physs(3, topo_id_internal) 
      max_phys(1) = ppm_max_physs(1, topo_id_internal)
      max_phys(2) = ppm_max_physs(2, topo_id_internal)
      max_phys(3) = ppm_max_physs(3, topo_id_internal)

#elif __KIND == __DOUBLE_PRECISION
      
      min_phys(1) = ppm_min_physd(1, topo_id_internal) 
      min_phys(2) = ppm_min_physd(2, topo_id_internal) 
      min_phys(3) = ppm_min_physd(3, topo_id_internal) 
      max_phys(1) = ppm_max_physd(1, topo_id_internal)
      max_phys(2) = ppm_max_physd(2, topo_id_internal)
      max_phys(3) = ppm_max_physd(3, topo_id_internal)
 
#endif


      IF (ppm_debug .GT. 0) THEN
         WRITE(mesg,'(A,3F15.3)' ) 'minimal extent', min_phys(1), min_phys(2), &
     &        min_phys(3)
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A,3F15.3)' ) 'maximal extent', max_phys(1), max_phys(2), &
     &        max_phys(3)
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
      ENDIF


      length_phys(1) = max_phys(1) - min_phys(1)
      length_phys(2) = max_phys(2) - min_phys(2)
      length_phys(3) = max_phys(3) - min_phys(3)


      !-------------------------------------------------------------------------
      ! Check if x_pencil topology
      !-------------------------------------------------------------------------
      Its_xpencil_topo = .TRUE.
      Its_xyslab_topo  = .TRUE.


      DO k=1,ppm_nsublist(topo_id_internal)
         idom = ppm_isublist(k,topo_id_internal)

#if   __KIND == __SINGLE_PRECISION
         length(1)=  ppm_max_subs(1,idom,topo_id_internal) - &
     &                                 ppm_min_subs(1,idom,topo_id_internal) 
         length(2)=  ppm_max_subs(2,idom,topo_id_internal) - &
     &                                 ppm_min_subs(2,idom,topo_id_internal) 
         IF( abs(length(1) - length_phys(1)).GT.(ppm_myepss) ) THEN 
            Its_xpencil_topo=.FALSE.
            Its_xyslab_topo =.FALSE.
         ENDIF
         IF( abs(length(2) - length_phys(2)).GT.(ppm_myepss) ) THEN 
            Its_xyslab_topo =.FALSE.
         ENDIF
#elif __KIND == __DOUBLE_PRECISION
         length(1)=ppm_max_subd(1,idom,topo_id_internal) - &
     &                                 ppm_min_subd(1,idom,topo_id_internal) 
         length(2)=ppm_max_subd(2,idom,topo_id_internal) - &
     &                                 ppm_min_subd(2,idom,topo_id_internal) 
         IF( abs(length(1) - length_phys(1)).GT.(ppm_myepsd) ) THEN 
            Its_xpencil_topo=.FALSE.
            Its_xyslab_topo =.FALSE.
         ENDIF
         IF( abs(length(2) - length_phys(2)).GT.(ppm_myepsd) ) THEN 
            Its_xyslab_topo =.FALSE.
         ENDIF
#endif

      ENDDO

      IF (ppm_debug .GT. 0) THEN
         IF(Its_xyslab_topo) THEN
            WRITE(mesg,'(A)' ) 'XY slab topology' 
            CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         
         ELSE
            WRITE(mesg,'(A)' ) 'Not XY slab topology' 
            CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         ENDIF
         IF(Its_xpencil_topo) THEN
            WRITE(mesg,'(A)' ) 'X pencil topology' 
            CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         
         ELSE
            WRITE(mesg,'(A)' ) 'Not X pencil topology' 
            CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         ENDIF
      ENDIF

#if __CASE == __SLAB

      !-------------------------------------------------------------------------
      !   Setting of xy slab topology  
      !-------------------------------------------------------------------------

      IF(Its_xyslab_topo) THEN

         topo_id_xpen = topo_id_user
         mesh_id_xpen = mesh_id_user
         topo_id_zpen = topo_ids(2)
         mesh_id_xpen_complex = mesh_ids(2)
         mesh_id_zpen      = mesh_ids(3)

      ELSE
         topo_id_xpen      = topo_ids(1)
         topo_id_zpen      = topo_ids(2)
         mesh_id_xpen         = mesh_ids(1)
         mesh_id_xpen_complex = mesh_ids(2)
         mesh_id_zpen         = mesh_ids(3)

      ENDIF

      IF (ppm_debug .GT. 0) THEN
         WRITE(mesg,'(A)' ) '  ID             topo  mesh' 
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A)' ) '-----------------------------------'
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'Original        ',topo_id_user, mesh_id_user  
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'XZ Slab        ', topo_id_xpen, mesh_id_xpen
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'XZ Slab  Complex', topo_id_xpen, &
     &       mesh_id_xpen_complex
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'Z Pencil Complex', topo_id_zpen, mesh_id_zpen
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
      ENDIF


      !-------------------------------------------------------------------------
      !  Decompose domain in xy slabs
      !-------------------------------------------------------------------------

      IF(.NOT.Its_xyslab_topo) THEN

         asign = ppm_param_assign_internal
         decomp = ppm_param_decomp_xy_slab

         CALL ppm_mktopo(xp,Npart,Nm,decomp,asign, min_phys,max_phys,  &
     &                    bcdef,ghostsize,topo_id_xpen,mesh_id_xpen,min_sub,&
     &                    max_sub,cost,sub2proc,nsubs,isublist,             &
     &                    nsublist,istart,ndata,info)
         
         topo_ids_tmp(1) = topo_id_user
         topo_ids_tmp(2) = topo_id_xpen

         mesh_ids_tmp(1) = mesh_id_user
         mesh_ids_tmp(2) = mesh_id_xpen


#if   __DIM == __SFIELD
         CALL ppm_fdsolver_map(DATA_fv,topo_ids_tmp, mesh_ids_tmp,          &
     &                          ghostsize, info)
#elif __DIM == __VFIELD
         CALL ppm_fdsolver_map(DATA_fv,lda_fv, topo_ids_tmp, mesh_ids_tmp,  &
     &                          ghostsize, info)
#endif


         topo_ids_tmp(1) = topo_id_xpen
         topo_ids_tmp(2) = topo_id_zpen

         topo_id_internal = ppm_field_topoid
         mesh_id_internal = ppm_meshid(topo_id_internal)%internal(mesh_id_xpen)
      ENDIF



#else

      !-------------------------------------------------------------------------
      !   Setting of x-pencil topology  
      !-------------------------------------------------------------------------

      IF(Its_xpencil_topo) THEN

         topo_id_xpen = topo_id_user
         mesh_id_xpen = mesh_id_user
         topo_id_ypen = topo_ids(2)
         topo_id_zpen = topo_ids(3)
         mesh_id_xpen_complex = mesh_ids(2)
         mesh_id_ypen      = mesh_ids(3)
         mesh_id_zpen      = mesh_ids(4)

      ELSE
         topo_id_xpen      = topo_ids(1)
         topo_id_ypen      = topo_ids(2)
         topo_id_zpen      = topo_ids(3)
         mesh_id_xpen         = mesh_ids(1)
         mesh_id_xpen_complex = mesh_ids(2)
         mesh_id_ypen         = mesh_ids(3)
         mesh_id_zpen         = mesh_ids(4)

      ENDIF


      IF (ppm_debug .GT. 0) THEN
         WRITE(mesg,'(A)' ) '  ID             topo  mesh' 
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A)' ) '-----------------------------------'
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'Original        ',topo_id_user, mesh_id_user  
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'X Pencil        ', topo_id_xpen, mesh_id_xpen
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'X Pencil Complex', topo_id_xpen, &
     &       mesh_id_xpen_complex
        CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
        WRITE(mesg,'(A,2I4)' ) 'Y Pencil Complex', topo_id_ypen, mesh_id_ypen
        CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
        WRITE(mesg,'(A,2I4)' ) 'Z Pencil Complex', topo_id_zpen, mesh_id_zpen
        CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
      ENDIF



      !-------------------------------------------------------------------------
      !  Decompose domain in xpencils
      !-------------------------------------------------------------------------

      IF(.NOT.Its_xpencil_topo) THEN

         asign = ppm_param_assign_internal
         decomp = ppm_param_decomp_xpencil

         CALL ppm_mktopo(xp,Npart,Nm,decomp,asign, min_phys,max_phys,  &
     &                    bcdef,ghostsize,topo_id_xpen,mesh_id_xpen,min_sub,&
     &                    max_sub,cost,sub2proc,nsubs,isublist,             &
     &                    nsublist,istart,ndata,info)


         topo_ids_tmp(1) = topo_id_user
         topo_ids_tmp(2) = topo_id_xpen

         mesh_ids_tmp(1) = mesh_id_user
         mesh_ids_tmp(2) = mesh_id_xpen


#if   __DIM == __SFIELD
         CALL ppm_fdsolver_map(DATA_fv,topo_ids_tmp, mesh_ids_tmp,         &
     &                              ghostsize, info)
#elif __DIM == __VFIELD
         CALL ppm_fdsolver_map(DATA_fv,lda_fv, topo_ids_tmp, mesh_ids_tmp,  & 
     &                             ghostsize, info)
#endif


         topo_ids_tmp(1) = topo_id_xpen
         topo_ids_tmp(2) = topo_id_ypen
         
         topo_id_internal = ppm_field_topoid
         mesh_id_internal = ppm_meshid(topo_id_internal)%internal(mesh_id_xpen)
      ENDIF

#endif

      !-------------------------------------------------------------------------
      ! Allocate complex array
      !-------------------------------------------------------------------------

      yhmax = 0
      zhmax = 0
      DO i=1,ppm_nsublist(topo_id_internal)
         idom = ppm_isublist(i,topo_id_internal)
         IF (ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes(2,idom)  &
    &           .GT.yhmax) THEN
           yhmax=ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes(2,idom)
         ENDIF
         IF (ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes(3,idom)  &
    &           .GT.zhmax) THEN
           zhmax=ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes(3,idom)
         ENDIF

      ENDDO

#if   __DIM == __SFIELD
      lda_DATA_fv_com(1)= Nm_com(1)
      lda_DATA_fv_com(2)= yhmax   
      lda_DATA_fv_com(3)= zhmax   
      lda_DATA_fv_com(4)= nsublist
#elif __DIM == __VFIELD
      lda_DATA_fv_com(1)= lda_fv
      lda_DATA_fv_com(2)= Nm_com(1)
      lda_DATA_fv_com(3)= yhmax   
      lda_DATA_fv_com(4)= zhmax   
      lda_DATA_fv_com(5)= nsublist
#endif

      iopt = ppm_param_alloc_fit
      CALL ppm_alloc(DATA_fv_com, lda_DATA_fv_com, iopt,info)


      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve',     &
     &        'data array',__LINE__,info)
         GOTO 9999
      ENDIF
     
 

      !-------------------------------------------------------------------------
      !  FFT - Transformation in x-direction
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      lda(1)=3
      lda(2)=ppm_nsubs(topo_id_internal)
      CALL ppm_alloc(ndata,lda,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve',     &
     &        'ndata array',__LINE__,info)
         GOTO 9999
      ENDIF
       

      ndata = ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes


      DO k=1,ppm_nsublist(topo_id_internal)

         idom = ppm_isublist(k,topo_id_internal)
         CALL  ppm_alloc(data_in, ndata(:,idom), iopt, info)
         
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve',     &
     &        'data_in array',__LINE__,info)
            GOTO 9999
         ENDIF


#if   __DIM == __SFIELD

         DO l=1, ndata(3,idom)
            DO j=1, ndata(2,idom)
               DO i=1, ndata(1,idom)
                  data_in(i,j,l) =  DATA_fv(i,j,l,k)
               ENDDO
            ENDDO
         ENDDO


         lda = ndata(:,idom)

#if __CASE == __SLAB
         CALL ppm_fdsolver_fft_fd_slab(data_in, lda,FFT_x,info) 
#elif __CASE == __INIT
         CALL ppm_fdsolver_fft_fd(data_in, lda,FFT_x,info) 
#else
         CALL ppm_util_fft_forward(data_in, lda,FFT_x,info) 
#endif

         DO l=1, lda(3)
            DO j=1, lda(2)
               DO i=1, lda(1)
                  DATA_fv_com(i,j,l,k) = FFT_x(i,j,l)
              ENDDO
           ENDDO
        ENDDO

         iopt = ppm_param_dealloc
         CALL  ppm_alloc(FFT_x, ndata(:,idom), iopt, info)
         IF (info .NE. 0) THEN
            WRITE(mesg,'(A)') 'could not deallocate memory'
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_solve',mesg,__LINE__,&
     &                                                                 info)
            GOTO 9999
         ENDIF


#elif __DIM == __VFIELD
         DO n=1,lda_fv  

            DO l=1, ndata(3,idom)
               DO j=1, ndata(2,idom)
                  DO i=1, ndata(1,idom)
                     data_in(i,j,l) =  DATA_fv(n,i,j,l,k)
                  ENDDO
               ENDDO
            ENDDO


            lda = ndata(:,idom)

#if __CASE == __SLAB
            CALL ppm_fdsolver_fft_fd_slab(data_in, lda,FFT_x,info) 
#elif __CASE == __INIT
            CALL ppm_fdsolver_fft_fd(data_in, lda,FFT_x,info) 
#else
            CALL ppm_util_fft_forward(data_in, lda,FFT_x,info) 
#endif
            DO l=1, lda(3)
               DO j=1, lda(2)
                  DO i=1, lda(1)
                     DATA_fv_com(n,i,j,l,k) = FFT_x(i,j,l)
                  ENDDO
               ENDDO
            ENDDO


         ENDDO

         iopt = ppm_param_dealloc
         CALL  ppm_alloc(data_in, ndata(:,idom), iopt, info)
         CALL  ppm_alloc(FFT_x, ndata(:,idom), iopt, info)


         IF (info .NE. 0) THEN
            WRITE(mesg,'(A)') 'could not deallocate memory'
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_solve',mesg,__LINE__,&
     &                                                                 info)
            GOTO 9999
         ENDIF

#endif



      ENDDO



      !-------------------------------------------------------------------------
      !  Decompose complex domain in xpencils
      !-------------------------------------------------------------------------
      CALL ppm_mesh_define(topo_id_xpen,Nm_com,min_phys,max_phys,             &
     &                      mesh_id_xpen_complex ,istart_xpen_complex,         &
     &                      ndata_xpen_complex, info)


      !-------------------------------------------------------------------------
      !  Decompose domain in ypencils
      !-------------------------------------------------------------------------
      decomp = ppm_param_decomp_ypencil
      CALL ppm_mktopo(xp,Npart,Nm_com,decomp,asign, min_phys,max_phys,  &
     &                bcdef,ghostsize,topo_id_ypen,mesh_id_ypen,min_sub,&
     &                max_sub,cost,sub2proc,nsubs,isublist,          &
     &                nsublist,istart_ypen,ndata_ypen,info)

      idom=isublist(1)
 
      !-------------------------------------------------------------------------
      !  Decompose domain in zpencils
      !-------------------------------------------------------------------------
      decomp = ppm_param_decomp_zpencil
      CALL ppm_mktopo(xp,Npart,Nm_com,decomp,asign, min_phys,max_phys,  &
     &                bcdef,ghostsize,topo_id_zpen,mesh_id_zpen,min_sub,&
     &                max_sub,cost,sub2proc,nsubs,isublist,          &
     &                nsublist,istart_zpen,ndata_zpen,info)
      idom = isublist(1)
 

#if __CASE == __SLAB
      GOTO 1000
#endif

      !-------------------------------------------------------------------------
      !  Transpose x-direction and y-direction
      !-------------------------------------------------------------------------
      topo_ids_tmp(1) = topo_id_xpen
      topo_ids_tmp(2) = topo_id_ypen

      mesh_ids_tmp(1) = mesh_id_xpen_complex
      mesh_ids_tmp(2) = mesh_id_ypen



#if   __DIM == __SFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,topo_ids_tmp, mesh_ids_tmp,        &
     &                                         ghostsize, info)
#elif __DIM == __VFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,lda_fv,topo_ids_tmp, mesh_ids_tmp, &
     &                ghostsize, info)
#endif
 

 

      DO k=1,ppm_nsublist(ppm_field_topoid)
 

         idom = ppm_isublist(k,ppm_field_topoid)     
         lda(1)=3
         lda(2)=idom
         iopt = ppm_param_alloc_fit
         CALL ppm_alloc(ndata_trans,lda,iopt, info)
         ndata_trans(1,idom)=ndata_ypen(2,idom)
         ndata_trans(2,idom)=ndata_ypen(1,idom)
         ndata_trans(3,idom)=ndata_ypen(3,idom)
         
#if  __DIM == __VFIELD
         DO n=1,lda_fv
#endif

            CALL ppm_alloc(data_com,ndata_trans(:,idom),iopt, info)

            DO l=1, ndata_ypen(3,idom)
               DO j=1, ndata_ypen(2,idom)
                  DO i=1, ndata_ypen(1,idom)

#if  __DIM == __SFIELD
                     data_com(j,i,l)= DATA_fv_com(i,j,l,k)
#elif __DIM == __VFIELD
                     data_com(j,i,l)= DATA_fv_com(n,i,j,l,k)
#endif
                  ENDDO
               ENDDO
            ENDDO
           
      !-------------------------------------------------------------------------
      !  FFT - Transformation in y-direction
      !-------------------------------------------------------------------------


#if __CASE == __INIT
            CALL ppm_fdsolver_fft_fd(data_com,ndata_trans(:,idom),FFT_xy,info) 
#else
            CALL ppm_util_fft_forward(data_com,ndata_trans(:,idom),FFT_xy,info) 
#endif

      !-------------------------------------------------------------------------
      !  Transpose back x-direction and y-direction
      !-------------------------------------------------------------------------
         DO l=1, ndata_ypen(3,idom)
            DO i=1, ndata_ypen(1,idom)
               DO j=1, ndata_ypen(2,idom)
#if  __DIM == __SFIELD
                   DATA_fv_com(i,j,l,k) = FFT_xy(j,i,l)
#elif __DIM == __VFIELD
                   DATA_fv_com(n,i,j,l,k) = FFT_xy(j,i,l)
#endif
                ENDDO
              ENDDO
           ENDDO

        ENDDO

        
        iopt = ppm_param_dealloc
        CALL ppm_alloc(data_com,ndata_trans(:,idom),iopt, info)
        CALL ppm_alloc(FFT_xy, ndata_trans(:,idom),iopt, info)
        IF (info .NE. 0) THEN
           WRITE(mesg,'(A)') 'could not deallocate memory'
           info = ppm_error_error
           CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_solve',mesg,__LINE__,&
     &                                                                 info)
           GOTO 9999
        ENDIF



#if  __DIM == __VFIELD
      ENDDO
#endif

 1000  CONTINUE! CASE SLAB continues



      !-------------------------------------------------------------------------
      !  Transpose x-direction and z-direction
      !-------------------------------------------------------------------------

      topo_ids_tmp(1) = topo_id_ypen 
      topo_ids_tmp(2) = topo_id_zpen 

      mesh_ids_tmp(1) = mesh_id_ypen
      mesh_ids_tmp(2) = mesh_id_zpen
 
#if __CASE == __SLAB
      topo_ids_tmp(1) = topo_id_xpen 
      mesh_ids_tmp(1) = mesh_id_xpen_complex
#endif


#if __DIM == __SFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,topo_ids_tmp, mesh_ids_tmp,        &
     &                                           ghostsize, info)
#elif __DIM == __VFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,lda_fv,topo_ids_tmp, mesh_ids_tmp, &
     & ghostsize, info)
#endif



      DO k=1,ppm_nsublist(ppm_field_topoid)
 
         idom = ppm_isublist(k,ppm_field_topoid)     
         lda(1)=3
         lda(2)=idom
         iopt = ppm_param_alloc_fit
         CALL ppm_alloc(ndata_trans,lda,iopt, info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve_3d',     &
   &                'ndata not allocated',__LINE__,info)
            GOTO 9999
         ENDIF
         ndata_trans(1,idom)=ndata_zpen(3,idom)
         ndata_trans(2,idom)=ndata_zpen(2,idom)
         ndata_trans(3,idom)=ndata_zpen(1,idom)

         CALL ppm_alloc(data_com,ndata_trans(:,idom),iopt, info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve_3d',     &
     &        'data_com not allocated',__LINE__,info)
            GOTO 9999
         ENDIF



#if __DIM == __VFIELD
         DO n=1,lda_fv
#endif

            DO l=1, ndata_zpen(3,idom)
               DO j=1, ndata_zpen(2,idom)
                  DO i=1, ndata_zpen(1,idom)

#if __DIM == __SFIELD
                     data_com(l,j,i)= DATA_fv_com(i,j,l,k)
#elif __DIM == __VFIELD
                     data_com(l,j,i)= DATA_fv_com(n,i,j,l,k)
#endif

                  ENDDO
               ENDDO
            ENDDO

          
      !-------------------------------------------------------------------------
      !  FFT - Transformation in z-direction
      !-------------------------------------------------------------------------
#if __CASE == __SLAB
           CALL ppm_fdsolver_fft_fd_z(data_com,ndata_trans(:,idom),FFT_xyz,info) 
#elif __CASE == __INIT
           CALL ppm_fdsolver_fft_fd_z(data_com,ndata_trans(:,idom),FFT_xyz,info) 
#else
           CALL ppm_util_fft_forward(data_com,ndata_trans(:,idom),FFT_xyz,info) 
#endif

      !-------------------------------------------------------------------------
      !  Solve Poisson Equation
      !-------------------------------------------------------------------------

            lda(1)=3
            lda(2)=idom
            iopt = ppm_param_alloc_fit
            CALL ppm_alloc(istart,lda,iopt,info)
            IF (info .NE. 0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve_3d',     &
   &           'istart not allocated',__LINE__,info)
               GOTO 9999
            ENDIF
          ! transpose x-z of istart
            istart(1,idom)= istart_zpen(3,idom)
            istart(2,idom)= istart_zpen(2,idom)
            istart(3,idom)= istart_zpen(1,idom)
            length(1)     = length_phys(3)
            length(2)     = length_phys(2)
            length(3)     = length_phys(1)
            Nm_poisson(1) = Nm_com(3)-1 ! corrected by -1 for ppm convention
            Nm_poisson(2) = Nm_com(2)-1 ! corrected by -1 for ppm convention
            Nm_poisson(3) = Nm_com(1)
            CALL ppm_fdsolver_poisson(FFT_xyz, ndata_trans(1:3,idom),       &
    &                             istart(1:3,idom),length, Nm_poisson, info)
            iopt = ppm_param_dealloc
            CALL ppm_alloc(istart,lda,iopt,info)
            IF (info .NE. 0) THEN
               info = ppm_error_fatal
               CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve_3d',     &
    &                        'istart not allocated',__LINE__,info)
               GOTO 9999
            ENDIF

      

      !-------------------------------------------------------------------------
      !  FFT - Backward Transformation in z-direction
      !-------------------------------------------------------------------------
#if __CASE == __SLAB
           CALL ppm_fdsolver_fft_bd_z(FFT_xyz,ndata_trans(:,idom),data_com,info) 
#elif __CASE == __INIT
           CALL ppm_fdsolver_fft_bd_z(FFT_xyz,ndata_trans(:,idom),data_com,info) 
#else
           CALL ppm_util_fft_backward(FFT_xyz,ndata_trans(:,idom),data_com,info) 
#endif

            iopt = ppm_param_dealloc
            CALL ppm_alloc(FFT_xyz,lda,iopt,info)


      !-------------------------------------------------------------------------
      !  Transpose back z-direction and x-direction
      !-------------------------------------------------------------------------
            DO i=1, ndata_zpen(1,idom)
               DO j=1, ndata_zpen(2,idom)
                  DO l=1, ndata_zpen(3,idom)

#if   __DIM == __SFIELD
                     DATA_fv_com(i,j,l,k) = data_com(l,j,i)
#elif   __DIM == __VFIELD
                     DATA_fv_com(n,i,j,l,k) = data_com(l,j,i)
#endif
                  ENDDO
               ENDDO
            ENDDO

#if   __DIM == __VFIELD
         ENDDO
#endif

         iopt = ppm_param_dealloc
         CALL ppm_alloc(data_com, lda, iopt,info)

         IF (info .NE. 0) THEN
            WRITE(mesg,'(A)') 'could not deallocate memory'
            info = ppm_error_error
            CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_solve',mesg,__LINE__,&
   &                                                                 info)
            GOTO 9999
         ENDIF


      ENDDO ! end of do-loop over k=1,nsublist

 
#if __CASE ==  __SLAB
      GOTO 2000
#endif

      topo_ids_tmp(1) = topo_id_zpen
      topo_ids_tmp(2) = topo_id_ypen
      
      mesh_ids_tmp(1) = mesh_id_zpen
      mesh_ids_tmp(2) = mesh_id_ypen


#if   __DIM == __SFIELD
      CALL ppm_fdsolver_map(DATA_fv_com, topo_ids_tmp, mesh_ids_tmp,       &
     &                     ghostsize, info)
#elif   __DIM == __VFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,lda_fv,topo_ids_tmp, mesh_ids_tmp, &
     &                     ghostsize, info)
#endif



      DO k=1,ppm_nsublist(ppm_field_topoid)
         idom = ppm_isublist(k,ppm_field_topoid)     
         iopt = ppm_param_alloc_fit
         ndata_trans(1,idom) = ndata_ypen(2,idom)
         ndata_trans(2,idom) = ndata_ypen(1,idom)
         ndata_trans(3,idom) = ndata_ypen(3,idom)

         CALL ppm_alloc(FFT_xy, ndata_trans(:,idom), iopt, info)

         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve',     &
     &         'FFT_xy array',__LINE__,info)
            GOTO 9999
         ENDIF
     

#if __DIM == __VFIELD         
         DO n=1,lda_fv
#endif


            DO l=1, ndata_ypen(3,idom) 
               DO j=1, ndata_ypen(2,idom)
                  DO i=1, ndata_ypen(1,idom)
#if __DIM == __SFIELD         
                     FFT_xy(j,i,l)= DATA_fv_com(i,j,l,k)
#elif __DIM == __VFIELD         
                     FFT_xy(j,i,l)= DATA_fv_com(n,i,j,l,k)
#endif
                  ENDDO
               ENDDO
            ENDDO


      !-------------------------------------------------------------------------
      !  FFT - Backward Transformation in y-direction
      !-------------------------------------------------------------------------

#if __CASE == __INIT
            CALL ppm_fdsolver_fft_bd(FFT_xy, ndata_trans(:,idom),data_com,info) 
#else
            CALL ppm_util_fft_backward(FFT_xy,ndata_trans(:,idom),data_com,info) 
#endif




      !-------------------------------------------------------------------------
      !  Transpose y-direction and x-direction
      !-------------------------------------------------------------------------

            DO l=1, ndata_ypen(3,idom)
               DO i=1, ndata_ypen(1,idom)
                  DO j=1, ndata_ypen(2,idom)

#if __DIM == __SFIELD         
                     DATA_fv_com(i,j,l,k) = data_com(j,i,l)
#elif __DIM == __VFIELD         
                     DATA_fv_com(n,i,j,l,k) = data_com(j,i,l)
#endif
                  ENDDO
               ENDDO
            ENDDO

#if __DIM == __VFIELD         
         ENDDO
#endif

         iopt = ppm_param_dealloc
         CALL ppm_alloc(data_com, lda, iopt,info)
         CALL ppm_alloc(FFT_xy, lda, iopt,info)
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve',     &
     &        'FFT_x array',__LINE__,info)
            GOTO 9999
         ENDIF


      ENDDO ! end of do-loop k=1,nsublist

 2000 CONTINUE ! CASE SLAB continues
      
      topo_ids_tmp(1) = topo_id_ypen
      topo_ids_tmp(2) = topo_id_xpen

      mesh_ids_tmp(1) = mesh_id_ypen
      mesh_ids_tmp(2) = mesh_id_xpen_complex


#if __CASE == __SLAB
      topo_ids_tmp(1) = topo_id_zpen
      mesh_ids_tmp(1) = mesh_id_zpen
#endif


#if __DIM == __SFIELD         
      CALL ppm_fdsolver_map(DATA_fv_com,topo_ids_tmp, mesh_ids_tmp,         &
     &                 ghostsize, info)     
#elif __DIM == __VFIELD         
      CALL ppm_fdsolver_map(DATA_fv_com,lda_fv,topo_ids_tmp, mesh_ids_tmp,  &
     &                   ghostsize, info)     
#endif
 






      DO k=1,ppm_nsublist(ppm_field_topoid)
         idom = ppm_isublist(k,topo_id_internal)

         iopt = ppm_param_alloc_fit
         CALL ppm_alloc(FFT_x,ndata_xpen_complex(:,idom), iopt, info)

 

#if __DIM == __VFIELD         
         DO n=1,lda_fv
#endif


            lda = ndata_xpen_complex(:,idom)

            DO l=1, ndata_xpen_complex(3,idom)
               DO j=1, ndata_xpen_complex(2,idom)
                  DO i=1, ndata_xpen_complex(1,idom)

#if __DIM == __SFIELD        
                     FFT_x(i,j,l) = DATA_fv_com(i,j,l,k)
#elif __DIM == __VFIELD         
                     FFT_x(i,j,l) = DATA_fv_com(n,i,j,l,k)
#endif

                  ENDDO
               ENDDO
            ENDDO

      !-------------------------------------------------------------------------
      !  FFT - Backward Transformation in x-direction
      !-------------------------------------------------------------------------

#if __CASE == __SLAB
            CALL ppm_fdsolver_fft_bd_slab(FFT_x, lda ,Result,info) 
#elif __CASE == __INIT
            CALL ppm_fdsolver_fft_bd(FFT_x, lda ,Result,info) 
#else
            CALL ppm_util_fft_backward(FFT_x, lda ,Result,info) 
#endif
          

      !-------------------------------------------------------------------------
      ! Correct Inverse by problem size factor 1/(Nx*Ny*Nz)
      ! Subtract 1 to fit ppm convention
      !-------------------------------------------------------------------------

            rN = 1.0_MK/real(((Nm(1)-1)*(Nm(2)-1)*(Nm(3)-1)), MK)

      
            DO l=1, lda(3)
               DO j=1, lda(2)
                  DO i=1, lda(1)

#if __DIM == __SFIELD     
                     DATA_fv(i,j,l,k) = Result(i,j,l)*rN
#elif __DIM == __VFIELD         
                     DATA_fv(n,i,j,l,k) = Result(i,j,l)*rN
#endif
                   
                  ENDDO
               ENDDO
            ENDDO


#if __DIM == __VFIELD         
         ENDDO
#endif

      ENDDO ! end of do-loop k=1,nsublist
 

      iopt = ppm_param_dealloc
      CALL ppm_alloc(FFT_x,lda,iopt, info)

#if __CASE == __SLAB
      !-------------------------------------------------------------------------
      ! Map to original topology if not xy-slab topology
      !-------------------------------------------------------------------------

      IF(.NOT.Its_xyslab_topo) THEN

         topo_ids_tmp(1) = topo_id_xpen
         topo_ids_tmp(2) = topo_id_user

         mesh_ids_tmp(1) = mesh_id_xpen
         mesh_ids_tmp(2) = mesh_id_user


#if __DIM == __SFIELD        
         CALL ppm_fdsolver_map(DATA_fv, topo_ids_tmp, mesh_ids_tmp,         &
    &                      ghostsize, info)
#elif __DIM == __VFIELD         
         CALL ppm_fdsolver_map(DATA_fv,lda_fv, topo_ids_tmp, mesh_ids_tmp,  &
    &                      ghostsize, info)
#endif


      ENDIF



#else 
      !-------------------------------------------------------------------------
      ! Map to original topology if not x-pencil topology
      !-------------------------------------------------------------------------

      IF(.NOT.Its_xpencil_topo) THEN

         topo_ids_tmp(1) = topo_id_xpen
         topo_ids_tmp(2) = topo_id_user

         mesh_ids_tmp(1) = mesh_id_xpen
         mesh_ids_tmp(2) = mesh_id_user


#if __DIM == __SFIELD        
         CALL ppm_fdsolver_map(DATA_fv, topo_ids_tmp, mesh_ids_tmp,        &
     &                                       ghostsize, info)
#elif __DIM == __VFIELD         
         CALL ppm_fdsolver_map(DATA_fv,lda_fv, topo_ids_tmp, mesh_ids_tmp, &
     &                                        ghostsize, info)
#endif


      ENDIF

#endif


      !-------------------------------------------------------------------------
      ! Deallocate memory
      !-------------------------------------------------------------------------


      iopt = ppm_param_dealloc
      CALL ppm_alloc(Result, lda, iopt,info)
      CALL ppm_alloc(data_in, lda, iopt,info)
      CALL ppm_alloc(data_com, lda, iopt,info)
      CALL ppm_alloc(FFT_x, lda, iopt,info)
      CALL ppm_alloc(FFT_xy, lda, iopt,info)
      CALL ppm_alloc(FFT_xyz, lda, iopt,info)
      CALL ppm_alloc(DATA_fv_com, lda_DATA_fv_com, iopt,info)
      CALL ppm_alloc(ndata,lda,iopt,info)
      CALL ppm_alloc(ndata_trans,lda,iopt, info)
      CALL ppm_alloc(min_sub,lda,iopt, info)
      CALL ppm_alloc(max_sub,lda,iopt, info)
      CALL ppm_alloc(cost,lda,iopt, info)
      CALL ppm_alloc(istart,lda,iopt, info)
      CALL ppm_alloc(istart_xpen_complex,lda,iopt, info)
      CALL ppm_alloc(istart_ypen,lda,iopt, info)
      CALL ppm_alloc(istart_zpen,lda,iopt, info)
      CALL ppm_alloc(istart_trans,lda,iopt, info)
      CALL ppm_alloc(ndata,lda,iopt, info)
      CALL ppm_alloc(ndata_xpen_complex,lda,iopt, info)
      CALL ppm_alloc(ndata_ypen,lda,iopt, info)
      CALL ppm_alloc(ndata_zpen,lda,iopt, info)
      CALL ppm_alloc(ndata_trans,lda,iopt, info)
      CALL ppm_alloc(sub2proc,lda,iopt, info)
      CALL ppm_alloc(isublist,lda,iopt, info)

      IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_solve',mesg,__LINE__,&
     &                                                                 info)
          GOTO 9999
      ENDIF



#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fdsolver_solve_3d',t0,info)
      RETURN



#if   __CASE == __SLAB
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_slab_3d_ss
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_slab_3d_sd
#endif
#endif



#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_slab_3d_vs
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_slab_3d_vd
#endif
#endif


#elif   __CASE == __INIT
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_init_3d_ss
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_init_3d_sd
#endif
#endif



#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_init_3d_vs
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_init_3d_vd
#endif
#endif


#else

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_3d_ss
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_3d_sd
#endif
#endif



#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_3d_vs
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_3d_vd
#endif
#endif
#endif
