#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_init
      !-------------------------------------------------------------------------
      !
      !  Purpose      : 
      !                initializes fieldsolver by creating all FFT-plans apriori
      !
      !  Input        :  
      !                  DATA_fv(:,:,:,:) (F) field data
      !                  lda_fv      (I) size of leading dimension in vector
      !                                  case
      !                  mesh_id_user(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) initial topology(xpencils)
      !                                  topo_ids(2) second  topology(ypencils)
      !                                  topo_ids(3) third   topology(zpencils)
      !                                         (3D only!!)
      !                                 
      !                  mesh_ids(3) (I) mesh IDs where the FFTs are performed
      !                                  mesh_ids(1) first mesh 
      !                                         (xpencils,real
      !                                  mesh_ids(2) second mesh 
      !                                         (xpencils,complex)
      !                                  mesh_ids(3) third mesh 
      !                                         (ypencils,complex)
      !                                  mesh_ids(4) forth mesh 
      !                                         (zpencils,complex) (3D only!!)
      !                  ghostsize(3) (I)ghostsize
      !
      !  Input/output : 
      !                                
      !                                
      !
      !  Output       :   
      !                 info          (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_init.f,v $
      !  Revision 1.13  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.11  2006/04/10 08:54:29  pchatela
      !  Made xp a REAL, DIMENSION(:,:), POINTER to get rid of warnings
      !
      !  Revision 1.10  2006/04/07 17:41:04  hiebers
      !  Changed type of variable xp to POINTER
      !
      !  Revision 1.9  2005/11/29 10:56:58  hiebers
      !  Changed arguments for fftwRoutines from INTEGER to INTEGER,DIMENSION(1)
      !  for  NAG compiler
      !
      !  Revision 1.8  2005/06/04 00:40:23  michaebe
      !  Cosmetics of cosmetics?
      !
      !  Revision 1.7  2005/06/04 00:37:49  michaebe
      !  cosmetics
      !
      !  Revision 1.6  2005/06/04 00:36:12  michaebe
      !  __ppm_module_fdsolver_init.f, line 491.25: 1516-023 (E) Subscript is
      !  out of bounds.
      !  This compiler warning eventually kinda got on my nerves so I acted.
      !
      !  Revision 1.5  2005/02/19 07:32:31  ivos
      !  Resolved CVS conflicts.
      !
      !  Revision 1.4  2005/02/18 08:01:55  hiebers
      !  minor changes in error messages
      !
      !  Revision 1.3  2005/02/16 22:22:59  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.2  2005/02/16 12:41:13  hiebers
      !  exchange FFTW_ESTIMATE by FFTW_MEASURE
      !
      !  Revision 1.1  2005/02/16 10:22:34  hiebers
      !  initial implementation
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_init_2d_sca_s(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_init_2d_sca_d(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_init_3d_sca_s(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_init_3d_sca_d(DATA_fv,mesh_id_user,topo_ids, &
     &  mesh_ids,ghostsize, info)
#endif
#endif
#endif

#if   __DIM == __VFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_init_2d_vec_s(DATA_fv,lda_fv,mesh_id_user, &
     &   topo_ids, mesh_ids,ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_init_2d_vec_d(DATA_fv,lda_fv,mesh_id_user, &
     &  topo_ids,mesh_ids,ghostsize,info)
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_init_3d_vec_s(DATA_fv,lda_fv,mesh_id_user, &
     &  topo_ids,mesh_ids,ghostsize,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_init_3d_vec_d(DATA_fv,lda_fv,mesh_id_user, &
     & topo_ids,mesh_ids,ghostsize, info)
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
      USE ppm_module_data
      USE ppm_module_data_mesh
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_mktopo
      USE ppm_module_data_fieldsolver
      USE ppm_module_error



      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION | __KIND == __COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif


      !-------------------------------------------------------------------------
      !  FFTW include
      !-------------------------------------------------------------------------
#ifdef  HAVE_LIBFFTW3
      INCLUDE "fftw3.f"
#endif
 

      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data
#if   __DIM == __SFIELD
#if   __MESH_DIM == __2D
      REAL(MK), DIMENSION(:,:,:),          POINTER   :: data_fv
#elif __MESH_DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:),        POINTER   :: data_fv
#endif
#endif


#if   __DIM == __VFIELD
#if   __MESH_DIM == __2D
      REAL(MK), DIMENSION(:,:,:,:),        POINTER   :: data_fv
      INTEGER,                             INTENT(IN):: lda_fv
#elif __MESH_DIM == __3D
      REAL(MK), DIMENSION(:,:,:,:,:),      POINTER   :: data_fv
      INTEGER,                             INTENT(IN):: lda_fv
#endif
#endif
      ! mesh ID of the data
      INTEGER                    , INTENT(IN)      :: mesh_id_user


      ! topo / mesh ID for the mapping
#if   __MESH_DIM == __2D
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: ghostsize
#elif __MESH_DIM == __3D
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(4)      , INTENT(IN   )   :: mesh_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: ghostsize
#endif

      INTEGER                    , INTENT(  OUT)   :: info


      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(ppm_kind_double)                   :: t0
      ! counters
      INTEGER                                 :: k, i, j
      CHARACTER(LEN=ppm_char)                 :: mesg
      

#if   __MESH_DIM == __2D
      REAL(MK), DIMENSION(:,:),              POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:),           POINTER   :: data_comp, data_compl
      INTEGER, DIMENSION(2)                            :: lda
#elif __MESH_DIM == __3D
      REAL(MK), DIMENSION(:,:,:),            POINTER   :: data_real
      COMPLEX(MK), DIMENSION(:,:,:),         POINTER   :: data_comp, data_compl
      INTEGER, DIMENSION(3)                            :: lda
#endif

      ! size of the data_in 
      INTEGER,DIMENSION(1)                    :: MB_in
      INTEGER                                 :: Nx_in, Ny_in, Nz_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out, Nz_out

#ifdef HAVE_LIBFFTW3
      INTEGER                            :: MBistride, MBrank, MBidist
      INTEGER, DIMENSION(1)              :: MBiembed, MBoembed
      INTEGER                            :: MBhowmany, MBodist
      INTEGER, DIMENSION(2)              :: iembed_slab, oembed_slab
#endif



#ifdef __MATHKEISAN
      ! MATHKEISAN variables for MathKeisan FFTs
      INTEGER                                 :: isign_fft, scale_fft
      INTEGER                                 :: incx, incy
 

      ! unused variables for initialization
      INTEGER                                 :: isys
      REAL(MK), DIMENSION(:),POINTER          :: work
#endif

      ! variables
      REAL(MK), DIMENSION(:,:),POINTER        :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp, asign
      REAL(MK), DIMENSION(3  )                :: min_phys, max_phys
      REAL(MK), DIMENSION(3  )                :: length
      REAL(MK), DIMENSION(3  )                :: length_phys
      INTEGER , DIMENSION(6  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id, mesh_id
      INTEGER                                 :: mesh_id_internal
      INTEGER                                 :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      INTEGER                                 :: mesh_id_zpen
      INTEGER                                 :: mesh_id_slab
      INTEGER                                 :: mesh_id_slab_complex
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart                       
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen
      INTEGER , DIMENSION(:,:), POINTER       :: istart_zpen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata 
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_zpen, ndata_slab
      INTEGER , DIMENSION(:  ), POINTER       :: sub2proc
      INTEGER , DIMENSION(:  ), POINTER       :: isublist
      INTEGER                                 :: nsublist, idom
      INTEGER                                 :: dim, n
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_user, topo_id_internal
      INTEGER                                 :: topo_id_xpen, topo_id_ypen
      INTEGER                                 :: topo_id_zpen
      INTEGER, DIMENSION(3)                   :: Nm, Nm_com
      INTEGER, DIMENSION(2)                   :: Nm_slab
      LOGICAL                                 :: Its_xpencil_topo
      LOGICAL                                 :: Its_xyslab_topo


      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_init',t0,info)

#if  !(defined(HAVE_LIBFFTW3) | defined(__MATHKEISAN))

      !-------------------------------------------------------------------------
      !  Error if FFTW library or NEC is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error


#ifndef HAVE_LIBFFTW3
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_init',  &
     &    'PPM was compiled without fftw support',__LINE__,info)
#endif

#ifndef __MATHKEISAN
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_init',  &
     &    'PPM was compiled without MATHKEISAN support',__LINE__,info)
#endif

      GOTO 9999   


#else

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
           IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_fdsolver_init',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (ppm_field_topoid .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_init',  &
     &            'No field topology has been defined so far',__LINE__,info)
              GOTO 9999
          ENDIF         
      ENDIF
      


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
             CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_init',      &
     &           mesg,__LINE__,info)
             GOTO 9999
         ENDIF
         lda(1) = ppm_dim
         lda(2) = ppm_nsubs(ppm_field_topoid)
         CALL ppm_alloc(ndata,lda,iopt,info)
         CALL ppm_alloc(ndata_slab,lda,iopt,info)
         IF (info .NE. 0) THEN
             info = ppm_error_fatal
             WRITE (mesg,'(A,I10,A)') 'allocating ',      &
     &                   ppm_nsublist(ppm_field_topoid) ,' ndata failed'
             CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_init',      &
     &           mesg,__LINE__,info)
             GOTO 9999
         ENDIF




      !-------------------------------------------------------------------------
      !  Initialize variables
      !-------------------------------------------------------------------------
      asign            = ppm_param_assign_internal
      topo_id_internal  = ppm_field_topoid
      topo_id_user      = ppm_user_topoid(topo_id_internal)
      mesh_id_internal  = ppm_meshid(topo_id_internal)%internal(mesh_id_user)
      nsublist          = ppm_nsublist(topo_id_internal)
      isublist          = ppm_isublist(:,topo_id_internal)
      Nm(1)             = ppm_cart_mesh(mesh_id_internal, topo_id_internal)%Nm(1)
      Nm(2)             = ppm_cart_mesh(mesh_id_internal, topo_id_internal)%Nm(2)
      Nm(3)             = ppm_cart_mesh(mesh_id_internal, topo_id_internal)%Nm(3)
      Nm_com(1)         = Nm(1)/2+1
      Nm_com(2)         = Nm(2)
      Nm_com(3)         = Nm(3)
      ndata             = ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes
      bcdef(1)          = ppm_param_bcdef_periodic 
      bcdef(2)          = ppm_param_bcdef_periodic 
      bcdef(3)          = ppm_param_bcdef_periodic 
      bcdef(4)          = ppm_param_bcdef_periodic 
      bcdef(5)          = ppm_param_bcdef_periodic 
      bcdef(6)          = ppm_param_bcdef_periodic 
!      

      Npart             = 0
      ! size of dummy arrays
      lda               = 1
      

#ifdef __MATHKEISAN
      CALL ppm_alloc(work,1,ppm_param_alloc_fit,info)
#endif
      CALL ppm_alloc(data_real,lda,ppm_param_alloc_fit,info)
      CALL ppm_alloc(data_comp,lda,ppm_param_alloc_fit,info)
      CALL ppm_alloc(data_compl,lda,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_init',     &
     &        'allocation failed for data_real, data_comp',__LINE__,info)
          GOTO 9999
      ENDIF

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
         CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
         WRITE(mesg,'(A,3F15.3)' ) 'maximal extent', max_phys(1), max_phys(2), &
     &        max_phys(3)
         CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
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
         IF(Its_xpencil_topo) THEN
            WRITE(mesg,'(A)' ) 'X pencil topology' 
            CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
         
         ELSE
            WRITE(mesg,'(A)' ) 'Not X pencil topology' 
            CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
         ENDIF
         IF(Its_xyslab_topo) THEN
            WRITE(mesg,'(A)' ) 'XY slab topology' 
            CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
         
         ELSE
            WRITE(mesg,'(A)' ) 'Not XY slab topology' 
            CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !   Setting of x-pencil topology  
      !-------------------------------------------------------------------------


      IF(Its_xpencil_topo) THEN

         topo_id_xpen = topo_id_user
         mesh_id_xpen = mesh_id_user
!         topo_ids(1)  = topo_id_xpen
         topo_id_ypen = topo_ids(2)
#if __MESH_DIM == __3D
         topo_id_zpen = topo_ids(3)
#endif
         mesh_id_xpen_complex = mesh_ids(2)
         mesh_id_ypen      = mesh_ids(3)
#if __MESH_DIM == __3D
         mesh_id_zpen      = mesh_ids(4)
#endif
      ELSE
         topo_id_xpen      = topo_ids(1)
         topo_id_ypen      = topo_ids(2)
#if __MESH_DIM == __3D
         topo_id_zpen      = topo_ids(3)
#endif         
         mesh_id_xpen         = mesh_ids(1)
         mesh_id_xpen_complex = mesh_ids(2)
         mesh_id_ypen         = mesh_ids(3)
#if __MESH_DIM == __3D
         mesh_id_zpen         = mesh_ids(4)
#endif
      ENDIF


      IF (ppm_debug .GT. 0) THEN
         WRITE(mesg,'(A)' ) '  ID             topo  mesh' 
         CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
         WRITE(mesg,'(A)' ) '-----------------------------------'
         CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'Original        ',topo_id_user, mesh_id_user  
         CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'X Pencil        ', topo_id_xpen, mesh_id_xpen
         CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'X Pencil Complex', topo_id_xpen, &
     &       mesh_id_xpen_complex
        CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
        WRITE(mesg,'(A,2I4)' ) 'Y Pencil Complex', topo_id_ypen, mesh_id_ypen
        CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
        WRITE(mesg,'(A,2I4)' ) 'Z Pencil Complex', topo_id_zpen, mesh_id_zpen
        CALL ppm_write(ppm_rank,'ppm_fdsolver_init',mesg,j)
      ENDIF



      !-------------------------------------------------------------------------
      !  Decompose domain in xy Slabs
      !-------------------------------------------------------------------------

      IF(.NOT.Its_xyslab_topo) THEN

         asign = ppm_param_assign_internal

         decomp = ppm_param_decomp_xy_slab
      
         CALL ppm_mktopo(xp,Npart,Nm,decomp,asign, min_phys,max_phys,  &
     &                    bcdef,ghostsize,topo_id_xpen,mesh_id_xpen,min_sub,&
     &                    max_sub,cost,sub2proc,nsubs,isublist,             &
     &                    nsublist,istart,ndata_slab,info)


      ELSE
         ndata_slab = ndata
      ENDIF




      !-------------------------------------------------------------------------
      !  Create Plan/ Table for xy slab topology
      !-------------------------------------------------------------------------
      

      idom = isublist(nsublist)
      ! substract 1 to fit ppm-convention

      Nx_in=ndata_slab(1,idom)-1
      Ny_in=ndata_slab(2,idom)-1
#if __MESH_DIM == __3D
      Nz_in=ndata_slab(3,idom)
#endif

      Nx_out = Nx_in/2 + 1
      Ny_out = Ny_in 
      Nz_out = Nz_in

      !-------------------------------------------------------------------------
      !  FFTW version xy slab
      !-------------------------------------------------------------------------


#ifdef HAVE_LIBFFTW3

#if __MESH_DIM == __3D
      MBRank    = 2
      MBIstride = 1

      MBHowmany = Nz_in
      lda(1) = Nx_in+1
      lda(2) = Ny_in+1
      lda(3) = Nz_in
 
      iopt   = ppm_param_alloc_fit

      CALL ppm_alloc(data_real,lda,iopt,info)
      data_real = 1.0_MK
      lda(1) = Nx_out
      lda(2) = Ny_out + 1
      lda(3) = Nz_out
      CALL ppm_alloc(data_comp,lda,iopt,info)
      data_comp = 1.0_MK

      Nm_slab(1) = Nx_in
      Nm_slab(2) = Ny_in


      iEmbed_slab(1)  = Nx_in+1
      iEmbed_slab(2)  = Ny_in+1
      oEmbed_slab(1)  = Nx_out 
      oEmbed_slab(2)  = Ny_out + 1
      MBiDist    = (Nx_in+1) *  (Ny_in+1)
      MBoDist    =  Nx_out * (Ny_out+1)
#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft_r2c(Plan_slab_fd_s, MBRank,Nm_slab(1), MBHowmany, &
           & data_real(1,1,1), iEmbed_slab(1),MBIstride,MBiDist, &
           & data_comp(1,1,1),oEmbed_slab(1),MBIstride,MBoDist,FFTW_MEASURE)
      CALL sfftw_execute(Plan_slab_fd_s)


#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft_r2c(Plan_slab_fd_d, MBRank,Nm_slab, MBHowmany, &
           & data_real(1,1,1), iEmbed_slab,MBIstride,MBiDist, &
           & data_comp(1,1,1), oEmbed_slab,MBIstride,MBoDist,FFTW_MEASURE)
#endif

      oEmbed_slab(1)  = Nx_in+1
      oEmbed_slab(2)  = Ny_in+1
      iEmbed_slab(1)  = Nx_out 
      iEmbed_slab(2)  = Ny_out + 1
      MBoDist    = (Nx_in+1) *  (Ny_in+1)
      MBiDist    =  Nx_out * (Ny_out+1)


#if   __KIND == __SINGLE_PRECISION

      CALL sfftw_plan_many_dft_c2r(Plan_slab_bd_s, MBRank,Nm_slab, MBHowmany, &
           & data_comp(1,1,1),iEmbed_slab,MBIstride,MBiDist, &
           & data_real(1,1,1),oEmbed_slab,MBIstride,MBoDist,FFTW_MEASURE)
      CALL sfftw_execute(Plan_slab_bd_s)

#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft_c2r(Plan_slab_bd_d, MBRank,Nm_slab, MBHowmany, &
           & data_comp(1,1,1), iEmbed_slab,MBIstride,MBiDist, &
           & data_real(1,1,1), oEmbed_slab,MBIstride,MBoDist,FFTW_MEASURE)
#endif


      CALL ppm_alloc(data_real,lda,ppm_param_dealloc,info)

#endif

#endif




      !-------------------------------------------------------------------------
      !  Decompose domain in xpencils
      !-------------------------------------------------------------------------

      IF(.NOT.Its_xpencil_topo) THEN

         asign = ppm_param_assign_internal

         ! xpencils decomposition
         decomp = ppm_param_decomp_xpencil

         CALL ppm_mktopo(xp,Npart,Nm,decomp,asign, min_phys,max_phys,  &
     &                    bcdef,ghostsize,topo_id_xpen,mesh_id_xpen,min_sub,&
     &                    max_sub,cost,sub2proc,nsubs,isublist,             &
     &                    nsublist,istart,ndata,info)


      ENDIF
      




      !-------------------------------------------------------------------------
      !  Create Plan/ Table for xpencil topology
      !-------------------------------------------------------------------------
      

      idom = isublist(nsublist)
      ! substract 1 to fit ppm-convention
      Nx_in=ndata(1,idom)-1

      Ny_in=ndata(2,idom)
#if __MESH_DIM == __3D
      Nz_in=ndata(3,idom)
#endif

      Nx_out = Nx_in/2 + 1
      Ny_out = Ny_in
      Nz_out = Nz_in

      !-------------------------------------------------------------------------
      !  FFTW version xpencil
      !-------------------------------------------------------------------------


#ifdef HAVE_LIBFFTW3

#if __MESH_DIM == __2D
      MBRank    = 1
      MBiEmbed(1)  = -1
      MBoEmbed(1)  = -1
      MBIstride = 1

      MBHowmany = Ny_in

      iopt   = ppm_param_alloc_fit

      lda(1) = Nx_in+1
      lda(2) = Ny_in
      CALL ppm_alloc(data_real,lda,iopt,info)
      lda(1) = Nx_out
      lda(2) = Ny_out
      CALL ppm_alloc(data_comp,lda,iopt,info)


#if   __KIND == __SINGLE_PRECISION
      MBiDist    = Nx_in+1
      MBoDist    = Nx_out
      MB_in (1)  = Nx_in
      CALL sfftw_plan_many_dft_r2c(Plan_fd_s, MBRank,MB_in, MBHowmany, &
           & data_real(1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_comp(1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
      MBiDist    = Nx_out
      MBoDist    = Nx_in+1

      CALL sfftw_plan_many_dft_c2r(Plan_bd_s, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_real(1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
#elif __KIND == __DOUBLE_PRECISION

      MBiDist    = Nx_in+1
      MBoDist    = Nx_out
      MB_in (1)  = Nx_in
      CALL dfftw_plan_many_dft_r2c(Plan_fd_d, MBRank,MB_in, MBHowmany, &
           & data_real(1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_comp(1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
      MBiDist    = Nx_out
      MBoDist    = Nx_in+1
      CALL dfftw_plan_many_dft_c2r(Plan_bd_d, MBRank,MB_in, MBHowmany, &
           & data_real(1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_comp(1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
#endif


#endif


#if __MESH_DIM == __3D
      MBRank    = 1
      MBiEmbed(1)  = -1
      MBoEmbed(1)  = -1
      MBIstride = 1

      MBHowmany = Ny_in*Nz_in

      iopt   = ppm_param_alloc_fit

      lda(1) = Nx_in+1
      lda(2) = Ny_in
      lda(3) = Nz_in
      CALL ppm_alloc(data_real,lda,iopt,info)

      lda(1) = Nx_out
      lda(2) = Ny_out
      lda(3) = Nz_out
      CALL ppm_alloc(data_comp,lda,iopt,info)

#if   __KIND == __SINGLE_PRECISION
      MBiDist    = Nx_in+1
      MBoDist    = Nx_out

      MB_in (1)  = Nx_in
      CALL sfftw_plan_many_dft_r2c(Plan_fd_s, MBRank,MB_in, MBHowmany, &
           & data_real(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_comp(1,1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
      MBiDist    = Nx_out
      MBoDist    = Nx_in+1

      CALL sfftw_plan_many_dft_c2r(Plan_bd_s, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_real(1,1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
#elif __KIND == __DOUBLE_PRECISION

      MBiDist    = Nx_in+1
      MBoDist    = Nx_out
      MB_in (1)  = Nx_in
      CALL dfftw_plan_many_dft_r2c(Plan_fd_d, MBRank,MB_in, MBHowmany, &
           & data_real(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_comp(1,1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
      MBiDist    = Nx_out
      MBoDist    = Nx_in+1
      CALL dfftw_plan_many_dft_c2r(Plan_bd_d, MBRank,MB_in, MBHowmany, &
           & data_real(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_comp(1,1,1),MBoEmbed(1),MBIstride,MBoDist,FFTW_MEASURE)
#endif


#endif

#endif


      !-------------------------------------------------------------------------
      ! MATHKEISAN version xpencil
      !-------------------------------------------------------------------------


#ifdef __MATHKEISAN

      lda_table = 2*Nx_in + 64
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(table_fd_s,lda_table,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_s,lda_table,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_init',     &
     &        'table_fd_s not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#elif __KIND == __DOUBLE_PRECISION
      CALL ppm_alloc(table_fd_d,lda_table,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_d,lda_table,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_init',     &
     &        'table_fd_d not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#endif

      scale_fft = 1
      isign_fft = 0


#if __MESH_DIM == __2D

#if   __KIND == __SINGLE_PRECISION
      CALL  scfft(isign_fft, Nx_in, scale_fft, data_real(1,1), &
     &                                 data_comp(1,1), table_fd_s, work, isys)
      CALL  csfft(isign_fft, Nx_in, scale_fft, data_comp(1,1), &
     &                                 data_real(1,1), table_bd_s, work, isys)

#elif __KIND == __DOUBLE_PRECISION
      CALL  dzfft(isign_fft, Nx_in, scale_fft, data_real(1,1), &
     &                                 data_comp(1,1), table_fd_d, work, isys)
      CALL  zdfft(isign_fft, Nx_in, scale_fft, data_comp(1,1), &
     &                                 data_real(1,1), table_bd_d, work, isys)
#endif

#endif




#if __MESH_DIM == __3D

#if   __KIND == __SINGLE_PRECISION
      CALL  scfft(isign_fft, Nx_in, scale_fft, data_real(1,1,1), &
     &                                 data_comp(1,1,1), table_fd_s, work, isys)
      CALL  csfft(isign_fft, Nx_in, scale_fft, data_comp(1,1,1), &
     &                                 data_real(1,1,1), table_bd_s, work, isys)

#elif __KIND == __DOUBLE_PRECISION
      CALL  dzfft(isign_fft, Nx_in, scale_fft, data_real(1,1,1), &
     &                                 data_comp(1,1,1), table_fd_d, work, isys)
      CALL  zdfft(isign_fft, Nx_in, scale_fft, data_comp(1,1,1), &
     &                                 data_real(1,1,1), table_bd_d, work, isys)
#endif

#endif



#endif







      !-------------------------------------------------------------------------
      !  Decompose domain in ypencils
      !-------------------------------------------------------------------------
      decomp = ppm_param_decomp_ypencil
      CALL ppm_mktopo(xp,Npart,Nm_com,decomp,asign, min_phys,max_phys,  &
     &                bcdef,ghostsize,topo_id_ypen,mesh_id_ypen,min_sub,&
     &                max_sub,cost,sub2proc,nsubs,isublist,          &
     &                nsublist,istart_ypen,ndata_ypen,info)



      !-------------------------------------------------------------------------
      !  Create Plan/ Table for ypencil topology
      !-------------------------------------------------------------------------


      idom = isublist(nsublist)
      ! substract 1 to fit ppm-convention
      Nx_in=ndata_ypen(2,idom)-1

      Ny_in=ndata_ypen(1,idom)

#if __MESH_DIM == __3D
      Nz_in=ndata_ypen(3,idom)
#endif

      Nx_out = Nx_in   
      Ny_out = Ny_in
      Nz_out = Nz_in

      !-------------------------------------------------------------------------
      !  FFTW version ypencil
      !-------------------------------------------------------------------------

#ifdef HAVE_LIBFFTW3

#if __MESH_DIM == __2D
      MBRank    = 1
      MBiEmbed(1)  = -1
      MBoEmbed(1)  = -1
      MBIstride = 1
      MBiDist    = Nx_in+1
      MBoDist    = Nx_out+1

      MB_in(1)  = Nx_in
      MBHowmany = Ny_in
      lda(1) = Nx_in+1
      lda(2) = Ny_in
      CALL ppm_alloc(data_comp,lda,iopt,info)
      CALL ppm_alloc(data_compl,lda,iopt,info)


#if   __KIND == __SINGLE_PRECISION


      CALL sfftw_plan_many_dft(Plan_fd_c_y, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_FORWARD, FFTW_MEASURE)
      CALL sfftw_plan_many_dft(Plan_bd_c_y, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_BACKWARD, FFTW_MEASURE)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft(Plan_fd_cc_y, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_FORWARD, FFTW_MEASURE)
      CALL dfftw_plan_many_dft(Plan_bd_cc_y, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_BACKWARD, FFTW_MEASURE)
#endif


#endif


#if __MESH_DIM == __3D
      MBRank    = 1
      MBiEmbed(1)  = -1
      MBoEmbed(1)  = -1
      MBIstride = 1
      MBiDist    = Nx_in+1
      MBoDist    = Nx_out+1

      MB_in(1)  = Nx_in
      MBHowmany = Ny_in*Nz_in
      lda(1) = Nx_in+1
      lda(2) = Ny_in
      lda(3) = Nz_in
      CALL ppm_alloc(data_comp,lda,iopt,info)
      CALL ppm_alloc(data_compl,lda,iopt,info)


#if   __KIND == __SINGLE_PRECISION


      CALL sfftw_plan_many_dft(Plan_fd_c_y, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_FORWARD, FFTW_MEASURE)
      CALL sfftw_plan_many_dft(Plan_bd_c_y, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_BACKWARD, FFTW_MEASURE)
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft(Plan_fd_cc_y, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_FORWARD, FFTW_MEASURE)
      CALL dfftw_plan_many_dft(Plan_bd_cc_y, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_BACKWARD, FFTW_MEASURE)
#endif

#endif

#endif


      !-------------------------------------------------------------------------
      ! MATHKEISAN version ypencil
      !-------------------------------------------------------------------------



#ifdef __MATHKEISAN

      lda_table_y = 2*Nx_in + 64
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(table_fd_c_y,lda_table_y,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_c_y,lda_table_y,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_init',     &
     &        'table_fd_s not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#elif __KIND == __DOUBLE_PRECISION
      CALL ppm_alloc(table_fd_cc_y,lda_table_y,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_cc_y,lda_table_y,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_init',     &
     &        'table_fd_d not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#endif

      scale_fft = 1
      isign_fft = 0
      incx      = 1
      incy      = 1



#if __MESH_DIM == __2D

#if   __KIND == __SINGLE_PRECISION
      CALL  cfft(isign_fft, Nx_in, scale_fft, data_comp(1,1), incx, &
     &            data_compl(1,1),incy,  table_fd_c_y, lda_table_y,work,1,isys)
#elif __KIND == __DOUBLE_PRECISION
      CALL  zfft(isign_fft, Nx_in, scale_fft, data_comp(1,1), incx, &
     &            data_compl(1,1),incy,  table_fd_cc_y, lda_table_y,work,1,isys)
#endif

#endif




#if __MESH_DIM == __3D

#if   __KIND == __SINGLE_PRECISION
      CALL  cfft(isign_fft, Nx_in, scale_fft, data_comp(1,1,1), incx, &
     &            data_compl(1,1,1),incy,  table_fd_c_y, lda_table_y,work,1,isys)
#elif __KIND == __DOUBLE_PRECISION
      CALL  zfft(isign_fft, Nx_in, scale_fft, data_comp(1,1,1), incx, &
     &            data_compl(1,1,1),incy,  table_fd_cc_y, lda_table_y,work,1,isys)
#endif

#endif



#endif





#if __MESH_DIM == __3D
 
      !-------------------------------------------------------------------------
      !  Decompose domain in zpencils
      !-------------------------------------------------------------------------

      decomp = ppm_param_decomp_zpencil
      CALL ppm_mktopo(xp,Npart,Nm_com,decomp,asign, min_phys,max_phys,  &
     &                bcdef,ghostsize,topo_id_zpen,mesh_id_zpen,min_sub,&
     &                max_sub,cost,sub2proc,nsubs,isublist,          &
     &                nsublist,istart_zpen,ndata_zpen,info)


      !-------------------------------------------------------------------------
      !  Create Plan/ Table for zpencil topology
      !-------------------------------------------------------------------------



      idom = isublist(nsublist)
      ! substract 1 to fit ppm-convention
      Nx_in=ndata_zpen(3,idom)-1
      Ny_in=ndata_zpen(2,idom)
      Nz_in=ndata_zpen(1,idom)

      Nx_out = Nx_in     
      Ny_out = Ny_in
      Nz_out = Nz_in

      !-------------------------------------------------------------------------
      !  FFTW version zpencil
      !-------------------------------------------------------------------------


#ifdef HAVE_LIBFFTW3


      MBRank    = 1
      MBiEmbed(1)  = -1
      MBoEmbed(1)  = -1
      MBIstride = 1
      MBiDist    = Nx_in+1
      MBoDist    = Nx_out+1

      MB_in(1)  = Nx_in
      MBHowmany = Ny_in*Nz_in
      lda(1) = Nx_in+1
      lda(2) = Ny_in
      lda(3) = Nz_in
      CALL ppm_alloc(data_comp,lda,iopt,info)
      CALL ppm_alloc(data_compl,lda,iopt,info)


#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_plan_many_dft(Plan_fd_c_z, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_FORWARD, FFTW_MEASURE)
      CALL sfftw_plan_many_dft(Plan_bd_c_z, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_BACKWARD, FFTW_MEASURE)

#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_plan_many_dft(Plan_fd_cc_z, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_FORWARD, FFTW_MEASURE)
      CALL dfftw_plan_many_dft(Plan_bd_cc_z, MBRank,MB_in, MBHowmany, &
           & data_comp(1,1,1), MBiEmbed(1),MBIstride,MBiDist, &
           & data_compl(1,1,1),MBoEmbed(1),MBIstride,MBoDist, &
           & FFTW_BACKWARD, FFTW_MEASURE)
#endif


#endif




      !-------------------------------------------------------------------------
      ! MATHKEISAN version zpencil
      !-------------------------------------------------------------------------




#ifdef __MATHKEISAN

      lda_table_z = 2*Nx_in + 64
#if   __KIND == __SINGLE_PRECISION
      CALL ppm_alloc(table_fd_c_z,lda_table_z,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_c_z,lda_table_z,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_init',     &
     &        'table_fd_s not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#elif __KIND == __DOUBLE_PRECISION
      CALL ppm_alloc(table_fd_cc_z,lda_table_z,ppm_param_alloc_fit,info)
      CALL ppm_alloc(table_bd_cc_z,lda_table_z,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_init',     &
     &        'table_fd_d not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
#endif

      scale_fft = 1
      isign_fft = 0
      incx      = 1
      incy      = 1

#if   __KIND == __SINGLE_PRECISION
      CALL  cfft(isign_fft, Nx_in, scale_fft, data_comp(1,1,1), incx, &
     &            data_compl(1,1,1),incy,  table_fd_c_z, lda_table_z,work,1,isys)
#elif __KIND == __DOUBLE_PRECISION
      CALL  zfft(isign_fft, Nx_in, scale_fft, data_comp(1,1), incx, &
     &            data_compl(1,1,1),incy,  table_fd_cc_z, lda_table_z,work,1,isys)
#endif



#endif




#endif

      !-------------------------------------------------------------------------
      !  Deallocate memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc


      CALL ppm_alloc(isublist,lda,iopt,info)
      CALL ppm_alloc(ndata,lda,iopt,info)
#ifdef __MATHKEISAN
      CALL ppm_alloc(work,1,iopt,info)
#endif
      CALL ppm_alloc(min_sub,lda,iopt,info)
      CALL ppm_alloc(max_sub,lda,iopt,info)
      CALL ppm_alloc(cost,lda,iopt,info)
      CALL ppm_alloc(istart,lda,iopt,info)
      CALL ppm_alloc(istart_ypen,lda,iopt,info)
      CALL ppm_alloc(istart_zpen,lda,iopt,info)
      CALL ppm_alloc(ndata,lda,iopt,info)
      CALL ppm_alloc(ndata_ypen,lda,iopt,info)
      CALL ppm_alloc(ndata_zpen,lda,iopt,info)
      CALL ppm_alloc(ndata_slab,lda,iopt,info)
      CALL ppm_alloc(sub2proc,lda,iopt,info)
      CALL ppm_alloc(isublist,lda,iopt,info)
      CALL ppm_alloc(data_real,lda,iopt,info)
      CALL ppm_alloc(data_comp,lda,iopt,info)
      CALL ppm_alloc(data_compl,lda,iopt,info)





#endif
     IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_init',mesg,__LINE__,&
     &                                                                 info)
          GOTO 9999
      ENDIF


 9999 CONTINUE
      CALL substop('ppm_fdsolver_init',t0,info)

      RETURN





#if   __DIM == __SFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_init_2d_sca_s 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_init_2d_sca_d 
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_init_3d_sca_s 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_init_3d_sca_d 
#endif
#endif
#endif

#if   __DIM == __VFIELD
#if   __MESH_DIM == __2D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_init_2d_vec_s 
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_init_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_init_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_init_3d_vec_d
#endif
#endif
#endif



