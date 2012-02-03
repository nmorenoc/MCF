#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_mg_alloc_bc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Does the (re)allocation of arrays of type
      !                 bc_value. It offers the same allocation
      !                 types as ppm_alloc for regular arrays.
      !
      !  Input        : lda(:)    (I) number of subdomains +
      !                               numbers of levels
      !                 iopt      (I) alloc action. One of:
      !                                 ppm_param_alloc_fit_preserve
      !                                 ppm_param_alloc_fit
      !                                 ppm_param_alloc_grow_preserve
      !                                 ppm_param_alloc_grow
      !                                 ppm_param_dealloc
      !
      !  Input/output : field (T) array of TYPE(bc_value)
      !                               which is to be (re)allocated.
      !
      !  Output       : info      (I) Return status. 0 if everything OK.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_mg_alloc_bc.f,v $
      !  Revision 1.7  2004/10/01 16:33:38  ivos
      !  cosmetics.
      !
      !  Revision 1.6  2004/10/01 16:09:10  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.5  2004/09/22 18:25:08  kotsalie
      !  MG new version
      !
      !  Revision 1.4  2004/07/26 13:49:18  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.3  2004/07/26 11:59:39  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.2  2004/07/26 07:46:37  ivos
      !  Changed to use single-interface modules. Updated all USE statements.
      !
      !  Revision 1.1  2004/06/29 14:43:14  kotsalie
      !  Needed for my type allocation
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if __DIM == __SFIELD
#if __MESH_DIM  == __2D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_alloc_bc_2d_sca_s(field,lda,iopt,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_alloc_bc_2d_sca_d(field,lda,iopt,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_alloc_bc_3d_sca_s(field,lda,iopt,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_alloc_bc_3d_sca_d(field,lda,iopt,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM  == __2D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_alloc_bc_2d_vec_s(field,lda,iopt,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_alloc_bc_2d_vec_d(field,lda,iopt,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_alloc_bc_3d_vec_s(field,lda,iopt,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_alloc_bc_3d_vec_d(field,lda,iopt,info)
#endif
#endif
#endif

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Module
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_mg
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error 
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND ==__SINGLE_PRECISION 
      TYPE(bc_value_2d_sca_s), DIMENSION(:), POINTER :: field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(bc_value_2d_sca_d), DIMENSION(:), POINTER :: field
#endif
#elif __MESH_DIM == __3D
#if __KIND ==__SINGLE_PRECISION
      TYPE(bc_value_3d_sca_s), DIMENSION(:), POINTER :: field 
#elif __KIND == __DOUBLE_PRECISION
      TYPE(bc_value_3d_sca_d), DIMENSION(:), POINTER  :: field
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND ==__SINGLE_PRECISION
      TYPE(bc_value_2d_vec_s), DIMENSION(:), POINTER :: field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(bc_value_2d_vec_d), DIMENSION(:), POINTER :: field
#endif
#elif __MESH_DIM == __3D
#if __KIND ==__SINGLE_PRECISION
      TYPE(bc_value_3d_vec_s), DIMENSION(:), POINTER :: field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(bc_value_3d_vec_d), DIMENSION(:), POINTER  :: field
#endif
#endif
#endif
      INTEGER                 , DIMENSION(:  ), INTENT(IN   ) :: lda
      INTEGER                                 , INTENT(IN   ) :: iopt

      INTEGER                           , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER            :: i
      INTEGER, DIMENSION(1) :: ldc
      REAL(MK)              :: t0

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND ==__SINGLE_PRECISION 
      TYPE(bc_value_2d_sca_s), DIMENSION(:), POINTER :: work_field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(bc_value_2d_sca_d), DIMENSION(:), POINTER :: work_field
#endif
#elif __MESH_DIM == __3D
#if __KIND ==__SINGLE_PRECISION
      TYPE(bc_value_3d_sca_s), DIMENSION(:), POINTER ::  work_field 
#elif __KIND == __DOUBLE_PRECISION
      TYPE(bc_value_3d_sca_d), DIMENSION(:), POINTER  :: work_field
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND ==__SINGLE_PRECISION
      TYPE(bc_value_2d_vec_s), DIMENSION(:), POINTER :: work_field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(bc_value_2d_vec_d), DIMENSION(:), POINTER :: work_field
#endif
#elif __MESH_DIM == __3D
#if __KIND ==__SINGLE_PRECISION
      TYPE(bc_value_3d_vec_s), DIMENSION(:), POINTER ::  work_field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(bc_value_3d_vec_d), DIMENSION(:), POINTER  :: work_field
#endif
#endif
#endif

      LOGICAL            :: lcopy,lalloc,lrealloc,ldealloc
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_mg_alloc_bc',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0 .AND. iopt .NE. ppm_param_dealloc) THEN
          IF (SIZE(lda,1) .LT. 1) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_alloc_field',  &
     &            'lda must be at least of length 1',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(1) .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_alloc_field',  &
     &            'lda(1) must be >= 0',__LINE__,info)
              GOTO 9999
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Check allocation type
      !-------------------------------------------------------------------------
      lcopy    = .FALSE.
      lalloc   = .FALSE.
      lrealloc = .FALSE.
      ldealloc = .FALSE.
      IF (iopt .EQ. ppm_param_alloc_fit_preserve) THEN
          !---------------------------------------------------------------------
          !  Fit memory and preserve the present contents
          !---------------------------------------------------------------------
          IF (ASSOCIATED(field)) THEN
              ldc(1) = SIZE(field,1)
              IF (ldc(1) .NE. lda(1))  THEN
                  lalloc   = .TRUE.
                  lrealloc = .TRUE.
                  lcopy    = .TRUE.
              ENDIF
          ELSE
              lalloc = .TRUE.
          ENDIF
      ELSEIF (iopt .EQ. ppm_param_alloc_fit) THEN
          !---------------------------------------------------------------------
          !  Fit memory and discard the present contents
          !---------------------------------------------------------------------
          IF (ASSOCIATED(field)) THEN
              ldc(1) = SIZE(field,1)
              IF (ldc(1) .NE. lda(1)) THEN
                  lalloc   = .TRUE.
                  lrealloc = .TRUE.
                  ldealloc = .TRUE.
              ENDIF
          ELSE
              lalloc = .TRUE.
          ENDIF
      ELSEIF (iopt .EQ. ppm_param_alloc_grow_preserve) THEN
          !---------------------------------------------------------------------
          !  Fit memory and preserve the present contents
          !---------------------------------------------------------------------
          IF (ASSOCIATED(field)) THEN
              ldc(1) = SIZE(field,1)
              IF (ldc(1) .LT. lda(1)) THEN
                  lalloc   = .TRUE.
                  lrealloc = .TRUE.
                  lcopy    = .TRUE.
              ENDIF
          ELSE
              lalloc = .TRUE.
          ENDIF
      ELSEIF (iopt .EQ. ppm_param_alloc_grow) THEN
          !---------------------------------------------------------------------
          !  Fit memory and discard the present contents
          !---------------------------------------------------------------------
          IF (ASSOCIATED(field)) THEN
              ldc(1) = SIZE(field,1)
              IF (ldc(1) .LT. lda(1)) THEN
                  lalloc   = .TRUE.
                  lrealloc = .TRUE.
                  ldealloc = .TRUE.
              ENDIF
          ELSE
              lalloc = .TRUE.
          ENDIF
      ELSEIF (iopt .EQ. ppm_param_dealloc) THEN
          !---------------------------------------------------------------------
          !  Deallocate
          !---------------------------------------------------------------------
          IF (ASSOCIATED(field)) THEN
              ldc(1) = SIZE(field,1)
              lrealloc = .TRUE.
              ldealloc = .TRUE.
          ENDIF
      ENDIF
              
      !-------------------------------------------------------------------------
      !  Perform the actual alloc action
      !-------------------------------------------------------------------------
      IF (lalloc) THEN
          !---------------------------------------------------------------------
          !  Allocate new array with new size and nullify all members
          !---------------------------------------------------------------------
          ALLOCATE(work_field(lda(1)),STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_mg_alloc_field',  &
     &            'new mgfield WORK_MESH',__LINE__,info)
              GOTO 9999
          ENDIF
              DO i=1,lda(1)
                  NULLIFY(work_field(i)%pbcvalue)
              ENDDO
      ENDIF

      IF (lcopy) THEN
          !---------------------------------------------------------------------
          !  Save the old contents
          !---------------------------------------------------------------------
              DO i=1,MIN(ldc(1),lda(1))
                  work_field(i)%pbcvalue => field(i)%pbcvalue
              ENDDO
      ENDIF

      IF (ldealloc) THEN
          !---------------------------------------------------------------------
          !  Deallocate the old contents
          !---------------------------------------------------------------------
              DO i=1,ldc(1)
                  IF (ASSOCIATED(field(i)%pbcvalue)) THEN
                      DEALLOCATE(field(i)%pbcvalue,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_mg_alloc_bc',&
     &                       'field correction FIELD%PBCVALUE',__LINE__,info)
                      ENDIF
                      NULLIFY(field(i)%pbcvalue)
                  ENDIF
              ENDDO
      ENDIF

      IF (lrealloc) THEN
          !---------------------------------------------------------------------
          !  Deallocate old pointer array
          !---------------------------------------------------------------------
          DEALLOCATE(field,STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_mg_alloc_bc',  &
     &            'mgfield data FIELD',__LINE__,info)
          ENDIF
          NULLIFY(field)
      ENDIF

      IF (lalloc) THEN
          !---------------------------------------------------------------------
          !  Point result to new array
          !---------------------------------------------------------------------
          field => work_field
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mg_alloc_bc',t0,info)
      RETURN

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_bc_2d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_bc_2d_sca_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_bc_3d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_bc_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_bc_2d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_bc_2d_vec_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_bc_3d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_bc_3d_vec_d
#endif
#endif
#endif

