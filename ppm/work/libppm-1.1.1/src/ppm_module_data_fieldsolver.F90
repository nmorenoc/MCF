#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      ! Module         :            ppm_module_data_fieldsolver
      !-------------------------------------------------------------------------
      !
      ! Purpose       :  data module of fieldsolver mostly containing the FFT
      !                        plans
      !               
      !
      ! Remarks       :
      !
      ! References    : 
      !
      ! Revisions     :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_fieldsolver.f,v $
      !  Revision 1.3  2005/02/17 17:47:00  hiebers
      !  Reimplementation
      !
      !
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

MODULE ppm_module_data_fieldsolver

   USE ppm_module_data,ONLY:ppm_kind_single,ppm_kind_double
   PRIVATE :: ppm_kind_single,ppm_kind_double

      ! FFTW Plans
      INTEGER*8 Plan_fd_s,     Plan_fd_d 
      INTEGER*8 Plan_slab_fd_s,Plan_slab_fd_d 
      INTEGER*8 Plan_fd_c_y,   Plan_fd_cc_y
      INTEGER*8 Plan_fd_c_z,   Plan_fd_cc_z 
      INTEGER*8 Plan_bd_s,     Plan_bd_d
      INTEGER*8 Plan_slab_bd_s,Plan_slab_bd_d 
      INTEGER*8 Plan_bd_c_y,   Plan_bd_cc_y
      INTEGER*8 Plan_bd_c_z,   Plan_bd_cc_z 


      ! MATHKEISAN variables for MathKeisan FFTs
      ! working storage
      REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_fd_s,   table_bd_s
      REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_fd_d,   table_bd_d
      REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_fd_c_y, table_bd_c_y
      REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_fd_cc_y,table_bd_cc_y
      REAL(ppm_kind_single), DIMENSION(:),POINTER :: table_fd_c_z, table_bd_c_z
      REAL(ppm_kind_double), DIMENSION(:),POINTER :: table_fd_cc_z,table_bd_cc_z

      ! the size of the working storage
      INTEGER, DIMENSION(1)              :: lda_table, lda_table_y, lda_table_z


END MODULE ppm_module_data_fieldsolver
