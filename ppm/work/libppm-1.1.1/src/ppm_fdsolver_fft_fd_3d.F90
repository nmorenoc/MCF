#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_fdsolver_fft_fd_3d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs Fast Fourier Transform forth
      !                 using the precomputed plans in ppm_fdsolver_init
      !
      !  Input        : data_in(:,:,:)  (F) 3d data array to be transformed
      !                                
      !  Input/output : lda(:)          (I) size of data
      !
      !  Output       : data_out(:,:,:) (F) transformed data
      !                 info            (I) return status. =0 if no error.
      !
      !  Remarks      : 
      !                                                  
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_fft_fd_3d.f,v $
      !  Revision 1.5  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.4  2005/02/17 17:50:26  hiebers
      !  removed typo in ppm_module_data_fieldsolver
      !
      !  Revision 1.3  2005/02/16 22:22:59  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.2  2005/02/16 12:41:44  hiebers
      !  removed print statements
      !
      !  Revision 1.1  2005/02/16 11:53:24  hiebers
      !  initial implementation
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __CASE == __SLAB

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_fft_fd_slab_3ds(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_fft_fd_slab_3dd(data_in,lda,data_out,info)
#endif


#else

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_fft_fd_3ds(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_fft_fd_3dd(data_in,lda,data_out,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_fdsolver_fft_fd_3dc(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_fdsolver_fft_fd_3dcc(data_in,lda,data_out,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z
      SUBROUTINE ppm_fdsolver_fft_fd_z_3dc(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX_Z
      SUBROUTINE ppm_fdsolver_fft_fd_z_3dcc(data_in,lda,data_out,info)
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
      USE ppm_module_data_fieldsolver
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION | __KIND ==__SINGLE_PRECISION_COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND ==__SINGLE_PRECISION_COMPLEX_Z 
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND ==__DOUBLE_PRECISION_COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_double
#elif __KIND ==__DOUBLE_PRECISION_COMPLEX_Z 
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
      ! input data
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:,:)    , INTENT(IN   ) :: data_in
#elif __KIND == __SINGLE_PRECISION_COMPLEX| __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:) , INTENT(IN   ) :: data_in
#elif __KIND==__SINGLE_PRECISION_COMPLEX_Z|__KIND ==__DOUBLE_PRECISION_COMPLEX_Z
      COMPLEX(MK), DIMENSION(:,:,:) , INTENT(IN   ) :: data_in
#endif
      ! size of array
      INTEGER, DIMENSION(3)           , INTENT(INOUT) :: lda
      ! output data, fast fourier transformed
      COMPLEX(MK), DIMENSION(:,:,:)   , POINTER       :: data_out
      INTEGER                         , INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                :: t0
      ! counters
      INTEGER                                 :: i,j,k,iopt
      ! size of the data_in 
      INTEGER                                 :: Nx_in, Ny_in, Nz_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out, Nz_out

#ifdef HAVE_LIBFFTW3
      ! FFTW Plan
      INTEGER*8                          :: Plan
      INTEGER                            :: mbistride, mbrank, mbidist, mbiembed
      INTEGER                            :: mboembed, mbhowmany, mbodist
#endif

#ifdef __MATHKEISAN
      ! MATHKEISAN variables for MathKeisan FFTs
      INTEGER                                 :: isign_fft,isys
#if __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      INTEGER                                 :: incx, incy
#elif __KIND==__SINGLE_PRECISION_COMPLEX_Z|__KIND==__DOUBLE_PRECISION_COMPLEX_Z
      INTEGER                                 :: incx, incy
#endif
      !scale of the transformation
      REAL(MK)                                :: scale_fft
      ! working storage
      REAL(MK), DIMENSION(:),POINTER          :: table, work
      ! the size of the working storage
      INTEGER, DIMENSION(1)                   :: lda_table, lda_work
#endif



      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_fft_fd_3d',t0,info)


#if  !(defined(HAVE_LIBFFTW3) | defined(__MATHKEISAN))

      !-------------------------------------------------------------------------
      !  Error if FFT library support is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error


#ifndef HAVE_LIBFFTW3
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_fft_fd_3d',  &
     &    'PPM was compiled without fftw support',__LINE__,info)
#endif

#ifndef __MATHKEISAN
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_fft_fd_3d',  &
     &    'PPM was compiled without MATHKEISAN support',__LINE__,info)
#endif

      GOTO 9999      
#else

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (SIZE(lda,1) .LT. 3) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_fft_fd_3d',  &
     &            'lda must be at least of size 3',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(1) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_fft_fd_3d',  &
     &            'mesh size: lda(1) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(2) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_fft_fd_3d',  &
     &            'mesh size: lda(2) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(3) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_fft_fd_3d',  &
     &            'mesh size: lda(3) must be >0' ,__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      ! Nx_in=lda(1)

      ! subtract 1 to fit ppm-convention
      Nx_in=lda(1)-1

      Ny_in=lda(2)
      Nz_in=lda(3)

      !-------------------------------------------------------------------------
      !  Allocate result array
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION        | __KIND == __DOUBLE_PRECISION
      Nx_out = Nx_in/2 + 1
#elif __KIND == __SINGLE_PRECISION_COMPLEX| __KIND == __DOUBLE_PRECISION_COMPLEX
      Nx_out = Nx_in     
#elif __KIND==__SINGLE_PRECISION_COMPLEX_Z|__KIND ==__DOUBLE_PRECISION_COMPLEX_Z
      Nx_out = Nx_in     
#endif
      Ny_out = Ny_in
      Nz_out = Nz_in

#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      lda(1) = Nx_out
#elif __KIND == __SINGLE_PRECISION_COMPLEX| __KIND == __DOUBLE_PRECISION_COMPLEX
      lda(1) = Nx_out+1      ! to fit ppm-convention
#elif __KIND==__SINGLE_PRECISION_COMPLEX_Z|__KIND ==__DOUBLE_PRECISION_COMPLEX_Z
      lda(1) = Nx_out+1      ! to fit ppm-convention
#endif


      lda(2) = Ny_out
      lda(3) = Nz_out
      CALL ppm_alloc(data_out,lda,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_fft_fd_3d',     &
     &        'fft result DATA_OUT',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  FFT - Transform in x-direction
      !-------------------------------------------------------------------------





      !-------------------------------------------------------------------------
      !  NEC version - Use MathKeisan Library
      !-------------------------------------------------------------------------
#ifdef __MATHKEISAN


#if __CASE == __SLAB
      !-------------------------------------------------------------------------
      !  not implemented yet
      !-------------------------------------------------------------------------

      info = ppm_error_fatal
      CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_fft_fd_2d',     &
     &        'version not implemented',__LINE__,info)
      GOTO 9999


#else

      !-------------------------------------------------------------------------
      !  Allocate working storage
      !-------------------------------------------------------------------------

      lda_work = 4*Nx_in
      CALL ppm_alloc(work,lda_work,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_fft_fd_2d',     &
     &        'work not allocated',__LINE__,info)
          GOTO 9999
      ENDIF
      

      !-------------------------------------------------------------------------
      !  Forward FFT 
      !-------------------------------------------------------------------------
      scale_fft = 1
      isign_fft = -1

      DO k=1,Nz_in
         DO j=1,Ny_in


#if   __KIND == __SINGLE_PRECISION
            CALL  scfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), &
     &                                 data_out(1,j,k), table_fd_s, work, isys)
#elif __KIND == __DOUBLE_PRECISION
            CALL  dzfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), &
     &                                 data_out(1,j,k), table_fd_d, work, isys)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
            CALL  cfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), incx, &
     &               data_out(1,j,k), incy, table_fd_c_y, lda_table_y, & 
     &              work, lda_work(1),isys)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
            CALL  zfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), incx, &
     &               data_out(1,j,k), incy, table_fd_cc_y, lda_table_y, & 
     &               work, lda_work(1),isys)
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z
            CALL  cfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), incx, &
     &               data_out(1,j,k), incy, table_fd_c_z, lda_table_z, & 
     &              work, lda_work(1),isys)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX_Z
            CALL  zfft(isign_fft, Nx_in, scale_fft, data_in(1,j,k), incx, &
     &               data_out(1,j,k), incy, table_fd_cc_z, lda_table_z, & 
     &               work, lda_work(1),isys)
#endif

         ENDDO      
      ENDDO

      !-------------------------------------------------------------------------
      !  Deallocate Memory
      !-------------------------------------------------------------------------
      CALL ppm_alloc(work,lda_work,ppm_param_dealloc,info)
      IF (info .NE. 0) THEN
         WRITE(mesg,'(A)') 'could not deallocate memory'
         info = ppm_error_error
         CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_fft_fd_3d',mesg,__LINE__,&
     &                                                                 info)
         GOTO 9999
      ENDIF
#endif

#else

      !-------------------------------------------------------------------------
      !  FFTW version for LINUX,...
      !-------------------------------------------------------------------------

#if __CASE == __SLAB

#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_execute_dft_r2c(Plan_slab_fd_s,data_in(1,1,1),data_out(1,1,1) )
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_execute_dft_r2c(Plan_slab_fd_d,data_in(1,1,1),data_out(1,1,1) )
#endif

#else

#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_execute_dft_r2c(Plan_fd_s,data_in(1,1,1),data_out(1,1,1) )
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_execute_dft_r2c(Plan_fd_d,data_in(1,1,1),data_out(1,1,1) )
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      CALL sfftw_execute_dft(Plan_fd_c_y,data_in(1,1,1),data_out(1,1,1) )
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      CALL dfftw_execute_dft(Plan_fd_cc_y,data_in(1,1,1),data_out(1,1,1) )
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z
      CALL sfftw_execute_dft(Plan_fd_c_z,data_in(1,1,1),data_out(1,1,1) )
#elif __KIND == __DOUBLE_PRECISION_COMPLEX_Z
      CALL dfftw_execute_dft(Plan_fd_cc_z,data_in(1,1,1),data_out(1,1,1) )
#endif
#endif

#endif


#endif 


#if __KIND == __SINGLE_PRECISION_COMPLEX | __KIND == __DOUBLE_PRECISION_COMPLEX
      !-------------------------------------------------------------------------
      !  Copy margin to conform with PPM convention
      !-------------------------------------------------------------------------
      DO j=1,Ny_out
         DO k=1,Nz_out
            data_out(lda(1),j,k) = data_out(1,j,k)
         ENDDO
      ENDDO     
#endif

#if __KIND ==__SINGLE_PRECISION_COMPLEX_Z| __KIND ==__DOUBLE_PRECISION_COMPLEX_Z
      !-------------------------------------------------------------------------
      !  Copy margin to conform with PPM convention
      !-------------------------------------------------------------------------
      DO j=1,Ny_out
         DO k=1,Nz_out
            data_out(lda(1),j,k) = data_out(1,j,k)
         ENDDO
      ENDDO     
#endif


#if __CASE == __SLAB
      DO i=1,Nx_out
         DO k=1,Nz_out
            data_out(i,lda(2),k) = data_out(i,1,k)
         ENDDO
      ENDDO     
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fdsolver_fft_fd_3d',t0,info)
      RETURN


#if   __CASE == __SLAB


#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_fft_fd_slab_3ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_fft_fd_slab_3dd
#endif


#else

#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_fft_fd_3ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_fft_fd_3dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_fdsolver_fft_fd_3dc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_fdsolver_fft_fd_3dcc
#elif __KIND == __SINGLE_PRECISION_COMPLEX_Z
      END SUBROUTINE ppm_fdsolver_fft_fd_z_3dc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX_Z
      END SUBROUTINE ppm_fdsolver_fft_fd_z_3dcc
#endif

#endif
