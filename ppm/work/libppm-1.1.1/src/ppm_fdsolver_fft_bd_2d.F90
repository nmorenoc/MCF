#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_fdsolver_fft_bd_2d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine performs Fast Fourier Transform backward
      !                 using the precomputed plans in ppm_fdsolver_init
      !
      !  Input        : data_in(:,:)   (F) 2d data array to be transformed 
      !
      !  Input/output : lda(:)         (I) size of data array
      !
      !  Output       : data_out(:,:)  (F) transformed data
      !                 info           (I) return status. =0 if no error.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fdsolver_fft_bd_2d.f,v $
      !  Revision 1.4  2006/09/04 18:34:42  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.3  2005/02/18 08:01:28  hiebers
      !  minor changes in error messages
      !
      !  Revision 1.2  2005/02/16 22:22:58  ivos
      !  Bugfix: replaced non-existing ppm_module_data_fdsolver with
      !  ppm_module_data_fieldsolver.
      !
      !  Revision 1.1  2005/02/16 11:54:02  hiebers
      !  initial implementation
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_fft_bd_2ds(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_fft_bd_2dd(data_in,lda,data_out,info)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      SUBROUTINE ppm_fdsolver_fft_bd_2dc(data_in,lda,data_out,info)
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      SUBROUTINE ppm_fdsolver_fft_bd_2dcc(data_in,lda,data_out,info)
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
#elif __KIND == __DOUBLE_PRECISION | __KIND ==__DOUBLE_PRECISION_COMPLEX 
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
      COMPLEX(MK), DIMENSION(:,:)       , INTENT(IN   ) :: data_in
      ! size of array
      INTEGER, DIMENSION(:)         , INTENT(INOUT) :: lda
      ! output data, inverse fast fourier transformed
#if   __KIND == __SINGLE_PRECISION        | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:)      , POINTER       :: data_out
#elif __KIND == __SINGLE_PRECISION_COMPLEX| __KIND == __DOUBLE_PRECISION_COMPLEX
      COMPLEX(MK), DIMENSION(:,:)   , POINTER       :: data_out 
#endif
      INTEGER                       , INTENT(  OUT) :: info

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                :: t0
      ! counters
      INTEGER                                 :: i,j,iopt
      ! size of the data_in 
      INTEGER                                 :: Nx_in, Ny_in
      ! size of the data_out 
      INTEGER                                 :: Nx_out, Ny_out

#ifdef HAVE_LIBFFTW3
      ! FFTW Plan
      INTEGER*8                               :: Plan
#endif

#ifdef __MATHKEISAN
      ! MATHKEISAN variables for MathKeisan FFTs
      INTEGER                                 :: isign_fft,isys
#if   __KIND == __SINGLE_PRECISION        | __KIND == __DOUBLE_PRECISION

#elif __KIND == __SINGLE_PRECISION_COMPLEX| __KIND == __DOUBLE_PRECISION_COMPLEX
      ! parameters for cfft (default=1)
      INTEGER                                 :: incx, incy
#endif

      ! scale_fft of the transformation
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
      CALL substart('ppm_fdsolver_fft_bd_2d',t0,info)

#if  !(defined(HAVE_LIBFFTW3) | defined(__MATHKEISAN))


      !-------------------------------------------------------------------------
      !  Error if library support is not available
      !-------------------------------------------------------------------------
      info = ppm_error_error


#ifndef HAVE_LIBFFTW3
      CALL ppm_error(ppm_err_nofftw,'ppm_fdsolver_fft_bd_2d',  &
     &    'PPM was compiled without fftw support',__LINE__,info)
#endif

#ifndef __MATHKEISAN
      CALL ppm_error(ppm_err_noMathKeisan,'ppm_fdsolver_fft_bd_2d',  &
     &    'PPM was compiled without MATHKEISAN support',__LINE__,info)
#endif


      GOTO 9999      
#else

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (SIZE(lda,1) .LT. 2) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_fft_bd_2d',  &
     &            'lda must be at least of size 2',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(1) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_fft_bd_2d',  &
     &            'mesh size: lda(1) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(2) .LE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_fft_bd_2d',  &
     &            'mesh size: lda(2) must be >0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF
      Nx_in=lda(1)
      Ny_in=lda(2)

      !-------------------------------------------------------------------------
      !  Allocate result array
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      Nx_out = (Nx_in-1)*2
#elif __KIND == __SINGLE_PRECISION_COMPLEX |__KIND == __DOUBLE_PRECISION_COMPLEX
      Nx_out = Nx_in -1
#endif
      Ny_out=Ny_in
      lda(1)=Nx_out+1
      lda(2)=Ny_out
      CALL ppm_alloc(data_out,lda,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_fft_bd_2d',     &
     &        'fft result DATA_OUT',__LINE__,info)
          GOTO 9999
      ENDIF
     
      !-------------------------------------------------------------------------
      !  FFT - Transform in x-direction
      !-------------------------------------------------------------------------



      !-------------------------------------------------------------------------
      !  NEC version - Use MathKeisan Library 1.5
      !-------------------------------------------------------------------------




#ifdef __MATHKEISAN

      !-------------------------------------------------------------------------
      !  Allocate working storage
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION         | __KIND == __DOUBLE_PRECISION
      lda_work(1) = 4*Nx_out
#elif __KIND == __SINGLE_PRECISION_COMPLEX | __KIND== __DOUBLE_PRECISION_COMPLEX
      lda_work(1) = 6*Nx_out
#endif

      CALL ppm_alloc(work,lda_work,ppm_param_alloc_fit,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_fdsolver_fft_bd_2d',     &
     &        'work not allocated',__LINE__,info)
          GOTO 9999
      ENDIF



      !-------------------------------------------------------------------------
      !  Forward FFT 
      !-------------------------------------------------------------------------


      scale_fft = 1
      isign_fft = 1
      DO j=1,Ny_in

#if   __KIND == __SINGLE_PRECISION
          CALL  csfft(isign_fft, Nx_out, scale_fft, data_in(1,j), &
     &                                 data_out(1,j), table_bd_s, work, isys)
#elif __KIND == __DOUBLE_PRECISION

          CALL  zdfft(isign_fft, Nx_out, scale_fft, data_in(1,j), &
     &                                 data_out(1,j), table_bd_d, work, isys)
#elif __KIND == __SINGLE_PRECISION_COMPLEX
          CALL  cfft(isign_fft, Nx_out, scale_fft, data_in(1,j), incx,        &
     &              data_out(1,j), incy, table_fd_c_y, lda_table_y, work, &
     &              lda_work(1), isys)

#elif __KIND == __DOUBLE_PRECISION_COMPLEX
          CALL  zfft(isign_fft, Nx_out, scale_fft, data_in(1,j), incx,        &
     &              data_out(1,j), incy, table_fd_c_y, lda_table_y, work, &
     &              lda_work(1), isys)
#endif



      ENDDO


      !-------------------------------------------------------------------------
      !  Deallocate Memory
      !-------------------------------------------------------------------------
      CALL ppm_alloc(work,lda_work,ppm_param_dealloc,info)

      IF (info .NE. 0) THEN
          WRITE(mesg,'(A)') 'could not deallocate memory'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fdsolver_bd_2d',mesg,__LINE__,&
     &                                                                 info)
          GOTO 9999
      ENDIF


#else

      !-------------------------------------------------------------------------
      !  FFTW version for LINUX,...
      !-------------------------------------------------------------------------


#if   __KIND == __SINGLE_PRECISION
      CALL sfftw_execute_dft_c2r(Plan_bd_s, data_in(1,1), data_out(1,1) )
#elif __KIND == __DOUBLE_PRECISION
      CALL dfftw_execute_dft_c2r(Plan_bd_d, data_in(1,1), data_out(1,1) )
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      CALL sfftw_execute_dft(Plan_bd_c_y,   data_in(1,1), data_out(1,1) )
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      CALL dfftw_execute_dft(Plan_bd_cc_y,   data_in(1,1), data_out(1,1) )
#endif


#endif
#endif


      !-------------------------------------------------------------------------
      !  Copy margin to conform with PPM convention
      !-------------------------------------------------------------------------
      DO j=1,Ny_out
            data_out(lda(1),j) = data_out(1,j)
      ENDDO     


      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fdsolver_fft_bd_2d',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_fft_bd_2ds
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_fft_bd_2dd
#elif __KIND == __SINGLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_fdsolver_fft_bd_2dc
#elif __KIND == __DOUBLE_PRECISION_COMPLEX
      END SUBROUTINE ppm_fdsolver_fft_bd_2dcc
#endif
