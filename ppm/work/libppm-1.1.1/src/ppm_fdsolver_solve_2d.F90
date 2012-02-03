#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_solve_2d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Poisson solver  using FFTW
      !                 Solves the negative poisson equation in a 2-dimensional 
      !                 periodic domain
      !                            - Laplacian of Phi = omega
      !                 It takes the field quantity omega in DATA_fv which is 
      !                 assumed to be periodic outside the domain. The most 
      !                 efficient initial topology for the fieldsolver is 
      !                 a x-pencil topology. After performing a FFT on the  
      !                 x-pencils the data is mapped onto y-pencils where 
      !                 again a FFT is performed.
      !                 The version solve_init is based on recomputed FFT-plans
      !                 created in ppm_fdsolver_init and destroyed in
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
      !   
      !                    ppm_fdsolver_finalize(info)
      !   
      !   
      !                 Standalone Usage:
      !
      !                    ppm_fdsolver_solve(arguments)
      !
      !  Input        : 
      !                  lda_fv      (F) leading dimension (vector case only)
      !                  mesh_id_user(I) mesh ID of the current data field mesh
      !                  topo_ids(2) (I) topology IDs on which the FFTs are 
      !                                  performed
      !                                  topo_ids(1) initial topology
      !                                               solve:      xpencils  
      !                                               solve_init: xpencils  
      !                                  topo_ids(2) second  topology
      !                                               solve:      ypencils  
      !                                               solve_init: ypencils  
      !
      !                  mesh_ids(3) (I) mesh IDs where the FFTs are performed
      !                                  mesh_ids(1) first mesh 
      !                                             solve:      xpencils,real  
      !                                             solve_init: xpencils,real  
      !                                  mesh_ids(2) second mesh 
      !                                             solve:      xpencils,complex
      !                                             solve_init: xpencils,complex  
      !                                  mesh_ids(3) third mesh 
      !                                             solve:      ypencils,complex
      !                                             solve_init: ypencils,complex  
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
      !  $Log: ppm_fdsolver_solve_2d.f,v $
      !  Revision 1.15  2006/09/04 18:34:44  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.13  2006/04/10 08:54:29  pchatela
      !  Made xp a REAL, DIMENSION(:,:), POINTER to get rid of warnings
      !
      !  Revision 1.12  2006/04/07 17:41:53  hiebers
      !  Changed type of variable xp to POINTER
      !
      !  Revision 1.11  2005/08/03 14:35:12  ivos
      !  Shortened the subrountine names to meet the F90 31-character limit
      !  standard. pgf90 on gonzales had problems...
      !
      !  Revision 1.10  2005/03/21 15:41:21  hiebers
      !  fixed bug in assigning lda_DATA_fv_com
      !
      !  Revision 1.9  2005/02/16 11:58:43  hiebers
      !  Major addition: finalized scalar and vector version
      !        implemented version <solve_init> that used FFT plans created
      !        in ppm_fdsolver_init and destroyed in ppm_fdsolver_finalize
      !        whereas the version <solve> remains a standalone solver
      !
      !  Revision 1.8  2004/11/03 11:14:16  hiebers
      !  Exchanged __SXF90 by __MATHKEISAN
      !
      !  Revision 1.7  2004/10/01 16:08:58  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.6  2004/08/19 13:14:09  hiebers
      !  debugged scalar/vector version
      !
      !  Revision 1.5  2004/08/17 12:12:17  hiebers
      !  added scalar/vector versions
      !
      !  Revision 1.4  2004/07/28 14:18:47  hiebers
      !  fixed bug in routine call of ppm_fdsolver_poisson
      !
      !  Revision 1.3  2004/07/26 15:38:47  ivos
      !  Inserted missing USE statements to resolve undefined references
      !  at link stage.
      !
      !  Revision 1.2  2004/07/26 11:59:38  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.1  2004/07/26 08:52:28  hiebers
      !  Recommited, formerly ppm_module_fieldsolver
      !
      !  Revision 1.2  2004/07/26 07:28:00  hiebers
      !  optimized loops
      !
      !  Revision 1.1  2004/05/19 15:32:20  hiebers
      !  implementation from scratch
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __CASE == __INIT
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_init_2d_ss(DATA_fv, mesh_id_user,&
         topo_ids, mesh_ids, ghostsize,  info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_init_2d_sd(DATA_fv, mesh_id_user,& 
         topo_ids, mesh_ids,  ghostsize,info)
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_init_2d_vs(DATA_fv, lda_fv, mesh_id_user, &
     &  topo_ids, mesh_ids, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_init_2d_vd(DATA_fv, lda_fv, mesh_id_user, &
     &  topo_ids, mesh_ids, ghostsize, info)
#endif
#endif


#else

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_2d_ss(DATA_fv, mesh_id_user,&
         topo_ids, mesh_ids, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_2d_sd(DATA_fv, mesh_id_user,& 
         topo_ids, mesh_ids, ghostsize, info)
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_2d_vs(DATA_fv, lda_fv, mesh_id_user, &
     &  topo_ids, mesh_ids, ghostsize, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_solve_2d_vd(DATA_fv, lda_fv, mesh_id_user, &
     &  topo_ids, mesh_ids, ghostsize, info)
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
      USE   ppm_module_fdsolver_fft_fd
      USE   ppm_module_fdsolver_fft_bd
      USE   ppm_module_fdsolver_poisson
  
      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif



      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      ! mesh ID of the data
      INTEGER                    , INTENT(IN)      :: mesh_id_user

#if   __DIM == __SFIELD
      ! data
      REAL(MK), DIMENSION(:,:,:),  POINTER         :: DATA_fv
#elif   __DIM == __VFIELD
      ! data
      REAL(MK), DIMENSION(:,:,:,:),  POINTER       :: DATA_fv
      INTEGER                    , INTENT(IN)      :: lda_fv
#endif

      ! topo / mesh ID for the mapping
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: topo_ids
      INTEGER, DIMENSION(3)      , INTENT(IN   )   :: mesh_ids 
      INTEGER, DIMENSION(2)      , INTENT(IN   )   :: ghostsize
      INTEGER                    , INTENT(  OUT)   :: info
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                :: t0
      ! parameter for alloc
      INTEGER, DIMENSION(2)                   :: lda
      ! counters
      INTEGER                                 :: k,i,j,n
      ! 1/ number of gridpoints
      REAL(MK)                                :: rN
      ! result array

#if   __DIM == __SFIELD
      INTEGER, PARAMETER                       :: lda_fv = 1
#endif

#if   __DIM == __SFIELD
      COMPLEX(MK), DIMENSION(:,:,:), POINTER   :: DATA_fv_com
      INTEGER , DIMENSION(3  )                 :: lda_DATA_fv_com
#elif   __DIM == __VFIELD
      COMPLEX(MK), DIMENSION(:,:,:,:), POINTER :: DATA_fv_com
      INTEGER , DIMENSION(4  )                 :: lda_DATA_fv_com
#endif


      REAL(MK), DIMENSION(:, :),       POINTER :: data_in 
      COMPLEX(MK), DIMENSION(:,:),     POINTER :: data_com
      COMPLEX(MK), DIMENSION(:,:),     POINTER :: FFT_x, FFT_xy
      REAL(MK), DIMENSION(:,:),        POINTER :: Result


      ! variables
      REAL(MK), DIMENSION(:,:),POINTER        :: xp
      INTEGER                                 :: Npart
      INTEGER                                 :: decomp, asign
      REAL(MK), DIMENSION(2  )                :: min_phys, max_phys
      REAL(MK), DIMENSION(2  )                :: length
      REAL(MK), DIMENSION(2  )                :: length_phys
      INTEGER , DIMENSION(4  )                :: bcdef 
      INTEGER                                 :: nsubs,topo_id, mesh_id
      INTEGER                                 :: mesh_id_internal
      INTEGER                                 :: mesh_id_xpen, mesh_id_ypen
      INTEGER                                 :: mesh_id_xpen_complex
      REAL(MK), DIMENSION(:,:), POINTER       :: min_sub,max_sub
      REAL(MK), DIMENSION(:  ), POINTER       :: cost
      INTEGER , DIMENSION(:,:), POINTER       :: istart, istart_xpen_complex  
      INTEGER , DIMENSION(:,:), POINTER       :: istart_ypen, istart_trans
      INTEGER , DIMENSION(:,:), POINTER       :: ndata, ndata_xpen_complex
      INTEGER , DIMENSION(:,:), POINTER       :: ndata_ypen, ndata_trans
      INTEGER , DIMENSION(:  ), POINTER       :: sub2proc
      INTEGER , DIMENSION(:  ), POINTER       :: isublist
      INTEGER                                 :: nsublist, idom
      INTEGER                                 :: dim, yhmax
      INTEGER                                 :: iopt
      INTEGER                                 :: topo_id_user, topo_id_internal
      INTEGER                                 :: topo_id_xpen, topo_id_ypen
      INTEGER, DIMENSION(2)                   :: topo_ids_tmp
      INTEGER, DIMENSION(3)                   :: mesh_ids_tmp 
      INTEGER, DIMENSION(2)                   :: Nm, Nm_com, Nm_poisson
      LOGICAL                                 :: Its_xpencil_topo
      CHARACTER(LEN=ppm_char)                 :: mesg


      
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_solve_2d',t0,info)


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




      !-------------------------------------------------------------------------
      !  Check if FFTW-Library is available of if NEC Library
      !-------------------------------------------------------------------------
#if  !(defined(HAVE_LIBFFTW3) | defined(__MATHKEISAN))
      info = ppm_error_error


#ifndef HAVE_LIBFFTW3
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_solve_2d',  &
     &            'fdsolver needs FFTW Library  ' ,__LINE__,info)
#endif

#ifndef __MATHKEISAN
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_solve_2d',  &
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
         CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve_2d',      &
     &           mesg,__LINE__,info)
         GOTO 9999
      ENDIF


      !-------------------------------------------------------------------------
      !  Initialize variables
      !-------------------------------------------------------------------------
      dim              = 1
      bcdef(1)         = ppm_param_bcdef_freespace 
      bcdef(2)         = ppm_param_bcdef_freespace
      bcdef(3)         = ppm_param_bcdef_freespace 
      bcdef(4)         = ppm_param_bcdef_freespace
      asign           = ppm_param_assign_internal
      topo_id_internal = ppm_field_topoid
      topo_id_user     = ppm_user_topoid(topo_id_internal)
      mesh_id_internal = ppm_meshid(topo_id_internal)%internal(mesh_id_user)
      topo_ids_tmp     = topo_ids
      mesh_ids_tmp     = mesh_ids
      nsublist         = ppm_nsublist(topo_id_internal)
      isublist         = ppm_isublist(:,topo_id_internal)
      Nm(1)            = ppm_cart_mesh(mesh_id_internal, topo_id_internal)%Nm(1)
      Nm(2)            = ppm_cart_mesh(mesh_id_internal, topo_id_internal)%Nm(2)
      Nm_com(1)        = Nm(1)/2+1
      Nm_com(2)        = Nm(2)
      Npart            = 0

#if   __KIND == __SINGLE_PRECISION
      
      min_phys(1) = ppm_min_physs(1, topo_id_internal) 
      min_phys(2) = ppm_min_physs(2, topo_id_internal) 
      max_phys(1) = ppm_max_physs(1, topo_id_internal)
      max_phys(2) = ppm_max_physs(2, topo_id_internal)

#elif __KIND == __DOUBLE_PRECISION
      
      min_phys(1) = ppm_min_physd(1, topo_id_internal) 
      min_phys(2) = ppm_min_physd(2, topo_id_internal) 
      max_phys(1) = ppm_max_physd(1, topo_id_internal)
      max_phys(2) = ppm_max_physd(2, topo_id_internal)
 
#endif

      
      IF (ppm_debug .GT. 0) THEN

         WRITE(mesg,'(A,2F15.3)' ) 'minimal extent', min_phys(1), min_phys(2) 
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A,2F15.3)' ) 'maximal extent', max_phys(1), max_phys(2) 
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
      ENDIF

      length_phys(1) = max_phys(1) - min_phys(1)
      length_phys(2) = max_phys(2) - min_phys(2)

      !-------------------------------------------------------------------------
      ! Check if x_pencil topology
      !-------------------------------------------------------------------------
      Its_xpencil_topo = .TRUE.


      DO k=1,ppm_nsublist(topo_id_internal)
          idom = ppm_isublist(k,topo_id_internal)

#if   __KIND == __SINGLE_PRECISION
         length(1)=ppm_max_subs(1,idom,topo_id_internal) &
     &                                 - ppm_min_subs(1,idom,topo_id_internal) 
         IF( abs(length(1) - length_phys(1)).GT.(ppm_myepss) ) THEN 
            Its_xpencil_topo=.FALSE.
         ENDIF
#elif __KIND == __DOUBLE_PRECISION
         length(1)=ppm_max_subd(1,idom,topo_id_internal) &
     &                                 - ppm_min_subd(1,idom,topo_id_internal) 
         IF( abs(length(1) - length_phys(1)).GT.(ppm_myepsd) ) THEN 
            Its_xpencil_topo=.FALSE.
         ENDIF
#endif

      ENDDO

      IF (ppm_debug .GT. 0) THEN

         IF(Its_xpencil_topo) THEN
            WRITE(mesg,'(A)' ) 'X pencil topology' 
            CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         
         ELSE
            WRITE(mesg,'(A)' ) 'Not X pencil topology' 
            CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !   Setting of x-pencil topology 
      !-------------------------------------------------------------------------

      IF(Its_xpencil_topo) THEN

         topo_id_xpen = topo_id_user
         mesh_id_xpen = mesh_id_user
         topo_ids_tmp(1)  = topo_id_xpen
         topo_id_ypen = topo_ids_tmp(2)
         mesh_id_xpen_complex = mesh_ids_tmp(2)
         mesh_id_ypen      = mesh_ids_tmp(1)

      ELSE
         topo_id_xpen      = topo_ids_tmp(1)
         topo_id_ypen      = topo_ids_tmp(2)
         mesh_id_xpen         = mesh_ids_tmp(1)
         mesh_id_xpen_complex = mesh_ids_tmp(2)
         mesh_id_ypen         = mesh_ids_tmp(3)

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
         WRITE(mesg,'(A,2I4)' ) 'X Pencil Complex', topo_id_xpen,              &
     &            mesh_id_xpen_complex
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
         WRITE(mesg,'(A,2I4)' ) 'Y Pencil Complex', topo_id_ypen, mesh_id_ypen
         CALL ppm_write(ppm_rank,'ppm_fdsolver_solve',mesg,j)
      ENDIF




      !-------------------------------------------------------------------------
      !  Decompose domain in xpencils
      !-------------------------------------------------------------------------

      IF(.NOT.Its_xpencil_topo) THEN

         asign = ppm_param_assign_internal

         ! xpencils decomposition
         decomp = ppm_param_decomp_xpencil

         CALL ppm_mktopo(xp,Npart,Nm,decomp,asign, min_phys,max_phys,        &
     &                    bcdef,ghostsize,topo_id_xpen,mesh_id_xpen,min_sub,&
     &                    max_sub,cost,sub2proc,nsubs,isublist,       &
     &                    nsublist,istart,ndata,info)


         topo_ids_tmp(1) = topo_id_user
         topo_ids_tmp(2) = topo_id_xpen

         mesh_ids_tmp(1) = mesh_id_user
         mesh_ids_tmp(2) = mesh_id_xpen

#if   __DIM == __SFIELD
         CALL ppm_fdsolver_map(DATA_fv, topo_ids_tmp, mesh_ids_tmp, info)
#elif __DIM == __VFIELD
         CALL ppm_fdsolver_map(DATA_fv,lda_fv, topo_ids_tmp, mesh_ids_tmp, info)
#endif
         topo_ids_tmp(1) = topo_id_xpen
         topo_ids_tmp(2) = topo_id_ypen

         topo_id_internal = ppm_field_topoid
         mesh_id_internal = ppm_meshid(topo_id_internal)%internal(mesh_id_xpen)
      ENDIF





      !-------------------------------------------------------------------------
      ! Allocate complex array
      !-------------------------------------------------------------------------

      yhmax = 0
      DO i=1,ppm_nsublist(topo_id_internal)
         idom = ppm_isublist(i,topo_id_internal)
         IF (ppm_cart_mesh(mesh_id_internal,topo_id_internal)%nnodes(2,idom)  &
    &           .GT.yhmax) THEN
            yhmax = ppm_cart_mesh(mesh_id_internal,topo_id_internal)          &
    &                    %nnodes(2,idom)
         ENDIF
      ENDDO


#if   __DIM == __SFIELD
      lda_DATA_fv_com(1)= Nm_com(1)
      lda_DATA_fv_com(2)= yhmax   
      lda_DATA_fv_com(3)= nsublist
#elif __DIM == __VFIELD
      lda_DATA_fv_com(1)= lda_fv
      lda_DATA_fv_com(2)= Nm_com(1)
      lda_DATA_fv_com(3)= yhmax   
      lda_DATA_fv_com(4)= nsublist
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
      lda(1)=2
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
         CALL  ppm_alloc(data_in, ndata(:,idom), ppm_param_alloc_fit, info)
         
         IF (info .NE. 0) THEN
            info = ppm_error_fatal
            CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_solve',     &
     &        'data_in array',__LINE__,info)
            GOTO 9999
         ENDIF
  
#if  __DIM == __VFIELD       
         DO n=1,lda_fv
#endif

            DO i=1, ndata(1,idom)
               DO j=1, ndata(2,idom)

#if  __DIM == __SFIELD
                  data_in(i,j) =  DATA_fv(i,j,k)
#elif __DIM == __VFIELD
                  data_in(i,j) =  DATA_fv(n,i,j,k)
#endif
               ENDDO
            ENDDO
   
#if __CASE == __INIT
            CALL ppm_fdsolver_fft_fd( data_in, ndata(:,idom), FFT_x, info)
#else 
            CALL ppm_util_fft_forward( data_in, ndata(:,idom), FFT_x, info)
#endif


            iopt = ppm_param_dealloc
            CALL  ppm_alloc(data_in, ndata(:,idom), iopt, info)


            DO i=1, ndata(1,idom)
               DO j=1, ndata(2,idom)
#if  __DIM == __SFIELD
                  DATA_fv_com(i,j,k) = FFT_x(i,j)
#elif __DIM == __VFIELD
                  DATA_fv_com(n,i,j,k) = FFT_x(i,j)
#endif
               ENDDO
            ENDDO

         ENDDO
   
#if __DIM == __VFIELD
      ENDDO
#endif     

      !-------------------------------------------------------------------------
      !  Decompose complex domain in xpencils
      !-------------------------------------------------------------------------
      CALL ppm_mesh_define(topo_id_xpen,Nm_com,min_phys,max_phys,      &
                          mesh_id_xpen_complex ,istart_xpen_complex,   &
     &                                      ndata_xpen_complex, info)


      !-------------------------------------------------------------------------
      !  Decompose domain in ypencils
      !-------------------------------------------------------------------------
      decomp = ppm_param_decomp_ypencil
      CALL ppm_mktopo(xp,Npart,Nm_com,decomp,asign, min_phys,max_phys,  &
     &                bcdef,ghostsize,topo_id_ypen,mesh_id_ypen,min_sub,&
     &                max_sub,cost,sub2proc,nsubs,isublist,          &
     &                nsublist,istart_ypen,ndata,info)
 


      !-------------------------------------------------------------------------
      !  Transpose x-direction and y-direction
      !-------------------------------------------------------------------------
      mesh_ids_tmp(1) = mesh_id_xpen_complex
      mesh_ids_tmp(2) = mesh_id_ypen
 
#if __DIM == __SFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,topo_ids_tmp, mesh_ids_tmp, info)
#elif __DIM == __VFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,lda_fv, topo_ids_tmp, mesh_ids_tmp, info)
#endif

      DO k=1,ppm_nsublist(ppm_field_topoid)
         idom = ppm_isublist(k,ppm_field_topoid)     
         lda(1)=2
         lda(2)=idom
         iopt = ppm_param_alloc_fit
         CALL ppm_alloc(ndata_trans,lda,iopt, info)
          
         ndata_trans(1,idom)=ndata(2,idom)
         ndata_trans(2,idom)=ndata(1,idom)
          
         CALL ppm_alloc(data_com,ndata_trans(:,idom),iopt, info)

#if __DIM == __VFIELD
         DO n=1,lda_fv
#endif
 
            DO i=1, ndata(1,idom)
               DO j=1, ndata(2,idom)
#if __DIM == __SFIELD
                  data_com(j,i)= DATA_fv_com(i,j,k)
#elif __DIM == __VFIELD
                  data_com(j,i)= DATA_fv_com(n,i,j,k)
#endif
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

            iopt = ppm_param_dealloc
            CALL ppm_alloc(data_com,ndata_trans(:,idom),iopt, info)

      
      !-------------------------------------------------------------------------
      !  Solve Poisson Equation
      !-------------------------------------------------------------------------
            ! transpose istart and physical length
            lda(1)=2
            lda(2)=idom
            iopt = ppm_param_alloc_fit
            CALL ppm_alloc(istart,lda,iopt,info)
            istart(1,idom)=istart_ypen(2,idom)
            istart(2,idom)=istart_ypen(1,idom)
            length(1)     = length_phys(2)
            length(2)     = length_phys(1)
            Nm_poisson(1) = Nm_com(2)-1 ! corrected by -1 for ppm convention
            Nm_poisson(2) = Nm_com(1)

            CALL ppm_fdsolver_poisson(FFT_xy, ndata_trans(1:2,idom),        &
     &                             istart(1:2,idom),length,Nm_poisson, info)
            iopt = ppm_param_dealloc
            CALL ppm_alloc(istart,lda,iopt,info)

      
      !-------------------------------------------------------------------------
      !  FFT - Backward Transformation in y-direction
      !-------------------------------------------------------------------------
#if __CASE == __INIT
            CALL ppm_fdsolver_fft_bd(FFT_xy, ndata_trans(:,idom), data_com,info) 
#else
            CALL ppm_util_fft_backward(FFT_xy, ndata_trans(:,idom), data_com,info) 
#endif

      !-------------------------------------------------------------------------
      !  Transpose y-direction and x-direction
      !-------------------------------------------------------------------------
            DO i=1, ndata(1,idom)
               DO j=1, ndata(2,idom)

#if __DIM == __SFIELD
                  DATA_fv_com(i,j,k) = data_com(j,i)
#elif __DIM == __VFIELD
                  DATA_fv_com(n,i,j,k) = data_com(j,i)
#endif
               ENDDO
            ENDDO

#if __DIM == __VFIELD
         ENDDO
#endif

         iopt = ppm_param_dealloc
         CALL ppm_alloc(data_com, lda, iopt,info)
      ENDDO ! end of do loop over k=1,nsublist
 

      topo_id    = topo_ids_tmp(1)
      topo_ids_tmp(1)= topo_ids_tmp(2)
      topo_ids_tmp(2)= topo_id
      mesh_id    = mesh_ids_tmp(1)
      mesh_ids_tmp(1)= mesh_ids_tmp(2)
      mesh_ids_tmp(2)= mesh_id

#if __DIM == __SFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,topo_ids_tmp, mesh_ids_tmp, info)
#elif __DIM == __VFIELD
      CALL ppm_fdsolver_map(DATA_fv_com,lda_fv,topo_ids_tmp, mesh_ids_tmp, info)
#endif

      DO k=1,ppm_nsublist(ppm_field_topoid)
         idom = ppm_isublist(k,ppm_field_topoid)

#if __DIM == __VFIELD
         DO n=1, lda_fv
#endif

            DO i=1, ndata_xpen_complex(1,idom)
               DO j=1, ndata_xpen_complex(2,idom)
#if __DIM == __SFIELD
                  FFT_x(i,j)= DATA_fv_com(i,j,k)
#elif __DIM == __VFIELD
                  FFT_x(i,j)= DATA_fv_com(n,i,j,k)
#endif
               ENDDO
            ENDDO
            

      !-------------------------------------------------------------------------
      !  FFT - Backward Transformation in y-direction
      !-------------------------------------------------------------------------
#if __CASE == __INIT
            CALL ppm_fdsolver_fft_bd(FFT_x,ndata_xpen_complex(:,idom),Result,info) 
#else
            CALL ppm_util_fft_backward(FFT_x,ndata_xpen_complex(:,idom),Result,info) 
#endif
            
      !-------------------------------------------------------------------------
      ! Correct Inverse by problem size factor 1/(Nx*Ny)     
      ! Subtract 1 to fit ppm convention
      !-------------------------------------------------------------------------

            rN = 1/dble((Nm(1)-1)*(Nm(2)-1))
            DO i=1,ndata_xpen_complex(1,idom)
               DO j=1,ndata_xpen_complex(2,idom)

#if __DIM == __SFIELD
                  DATA_fv(i,j,k)= Result(i,j)*rN
#elif __DIM == __VFIELD
                  DATA_fv(n,i,j,k)= Result(i,j)*rN
#endif
               ENDDO
            ENDDO
#if __DIM == __VFIELD
         ENDDO
#endif

      ENDDO ! end of do loop k=1,nsublist
      

      !-------------------------------------------------------------------------
      ! Map to original topology if not x-pencil topology
      !-------------------------------------------------------------------------

      IF(.NOT.Its_xpencil_topo) THEN
         topo_ids_tmp(1) = topo_ids_tmp(2)
         topo_ids_tmp(2) = topo_id_user
         mesh_ids_tmp(1) = mesh_id_xpen
         mesh_ids_tmp(2) = mesh_id_user

#if __DIM == __SFIELD
         CALL ppm_fdsolver_map(DATA_fv, topo_ids_tmp, mesh_ids_tmp, info)
#elif __DIM == __VFIELD
         CALL ppm_fdsolver_map(DATA_fv, lda_fv, topo_ids_tmp, mesh_ids_tmp, info)
#endif
      ENDIF 

      !-------------------------------------------------------------------------
      ! Deallocate memory
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(DATA_fv_com, lda_DATA_fv_com, iopt,info)
      CALL ppm_alloc(ndata,lda,iopt,info)
      CALL ppm_alloc(ndata_trans,lda,iopt, info)
      CALL ppm_alloc(data_in,lda,iopt, info)
      CALL ppm_alloc(FFT_x,lda,iopt, info)
      CALL ppm_alloc(FFT_xy,lda,iopt, info)
      CALL ppm_alloc(min_sub,lda,iopt, info)
      CALL ppm_alloc(max_sub,lda,iopt, info)
      CALL ppm_alloc(cost,lda,iopt, info)
      CALL ppm_alloc(istart,lda,iopt, info)
      CALL ppm_alloc(istart_xpen_complex,lda,iopt, info)
      CALL ppm_alloc(istart_ypen,lda,iopt, info)
      CALL ppm_alloc(istart_trans,lda,iopt, info)
      CALL ppm_alloc(ndata,lda,iopt, info)
      CALL ppm_alloc(ndata_xpen_complex,lda,iopt, info)
      CALL ppm_alloc(ndata_ypen,lda,iopt, info)
      CALL ppm_alloc(ndata_trans,lda,iopt, info)
      CALL ppm_alloc(sub2proc,lda,iopt, info)
      CALL ppm_alloc(isublist,lda,iopt, info)

      IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_solve2d',mesg,__LINE__,&
     &                                                                 info)
          GOTO 9999
      ENDIF


#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fdsolver_solve_2d',t0,info)
      RETURN

#if   __CASE == __INIT
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_init_2d_ss
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_init_2d_sd
#endif
#endif


#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_init_2d_vs
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_init_2d_vd
#endif
#endif

#else
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_2d_ss
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_2d_sd
#endif
#endif


#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_2d_vs
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_solve_2d_vd
#endif
#endif

#endif
