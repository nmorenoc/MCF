#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_mg_alloc_field
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Does the (re)allocation of arrays of type
      !                 mg_field. It offers the same allocation
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
      !  Input/output : field (T) array of TYPE(mg_field)
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
      !  $Log: ppm_mg_alloc_field.f,v $
      !  Revision 1.8  2006/07/21 11:30:57  kotsalie
      !  FRIDAY
      !
      !  Revision 1.7  2004/10/01 16:33:39  ivos
      !  cosmetics.
      !
      !  Revision 1.6  2004/10/01 16:09:10  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.5  2004/09/22 18:24:06  kotsalie
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
      !  Revision 1.1  2004/06/29 14:42:16  kotsalie
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
      SUBROUTINE ppm_mg_alloc_field_2d_sca_s(field,lda,iopt,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_alloc_field_2d_sca_d(field,lda,iopt,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_alloc_field_3d_sca_s(field,lda,iopt,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_alloc_field_3d_sca_d(field,lda,iopt,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM  == __2D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_alloc_field_2d_vec_s(field,lda,iopt,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_alloc_field_2d_vec_d(field,lda,iopt,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_mg_alloc_field_3d_vec_s(field,lda,iopt,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_mg_alloc_field_3d_vec_d(field,lda,iopt,info)
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
      INTEGER                 , DIMENSION(:  ), INTENT(IN   ) :: lda
      INTEGER                                 , INTENT(IN   ) :: iopt
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND ==__SINGLE_PRECISION 
      TYPE(mg_field_2d_sca_s), DIMENSION(:,:), POINTER :: field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_2d_sca_d), DIMENSION(:,:), POINTER :: field
#endif
#elif __MESH_DIM == __3D
#if __KIND ==__SINGLE_PRECISION
      TYPE(mg_field_3d_sca_s), DIMENSION(:,:), POINTER :: field 
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_3d_sca_d),DIMENSION(:,:), POINTER  :: field
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND ==__SINGLE_PRECISION
      TYPE(mg_field_2d_vec_s), DIMENSION(:,:), POINTER :: field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_2d_vec_d), DIMENSION(:,:), POINTER :: field
#endif
#elif __MESH_DIM == __3D
#if __KIND ==__SINGLE_PRECISION
      TYPE(mg_field_3d_vec_s), DIMENSION(:,:), POINTER :: field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_3d_vec_d),DIMENSION(:,:), POINTER  :: field
#endif
#endif
#endif


      INTEGER                           , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER            :: i,j
      INTEGER, DIMENSION(2) :: ldc
      REAL(MK)              :: t0

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND ==__SINGLE_PRECISION
      TYPE(mg_field_2d_sca_s), DIMENSION(:,:), POINTER :: work_field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_2d_sca_d), DIMENSION(:,:), POINTER :: work_field
#endif
#elif __MESH_DIM == __3D
#if __KIND ==__SINGLE_PRECISION
      TYPE(mg_field_3d_sca_s), DIMENSION(:,:), POINTER :: work_field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_3d_sca_d), DIMENSION(:,:), POINTER :: work_field
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND ==__SINGLE_PRECISION
      TYPE(mg_field_2d_vec_s), DIMENSION(:,:), POINTER :: work_field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_2d_vec_d), DIMENSION(:,:), POINTER :: work_field
#endif
#elif __MESH_DIM == __3D
#if __KIND ==__SINGLE_PRECISION
      TYPE(mg_field_3d_vec_s), DIMENSION(:,:), POINTER :: work_field
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_3d_vec_d), DIMENSION(:,:), POINTER :: work_field
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
      CALL substart('ppm_mg_alloc_field',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0 .AND. iopt .NE. ppm_param_dealloc) THEN
          IF (SIZE(lda,1) .LT. 2) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_alloc_field',  &
     &            'lda must be at least of length 2',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(1) .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_alloc_field',  &
     &            'lda(1) must be >= 0',__LINE__,info)
              GOTO 9999
          ENDIF
          IF (lda(2) .LT. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_mg_alloc_field',  &
     &            'lda(2) must be >= 0',__LINE__,info)
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
              ldc(2) = SIZE(field,2)
              IF ((ldc(1) .NE. lda(1)) .OR. (ldc(2) .NE. lda(2))) THEN
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
              ldc(2) = SIZE(field,2)
              IF ((ldc(1) .NE. lda(1)) .OR. (ldc(2) .NE. lda(2))) THEN
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
              ldc(2) = SIZE(field,2)
              IF ((ldc(1) .LT. lda(1)) .OR. (ldc(2) .LT. lda(2))) THEN
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
              ldc(2) = SIZE(field,2)
              IF ((ldc(1) .LT. lda(1)) .OR. (ldc(2) .LT. lda(2))) THEN
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
              ldc(2) = SIZE(field,2)
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
          ALLOCATE(work_field(lda(1),lda(2)),STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_mg_alloc_field',  &
     &            'new mgfield WORK_MESH',__LINE__,info)
              GOTO 9999
          ENDIF
          DO j=1,lda(2)
              DO i=1,lda(1)
                  NULLIFY(work_field(i,j)%uc)
                  NULLIFY(work_field(i,j)%fc)
                  NULLIFY(work_field(i,j)%err)
                  NULLIFY(work_field(i,j)%bcvalue)
              ENDDO
          ENDDO
      ENDIF

      IF (lcopy) THEN
          !---------------------------------------------------------------------
          !  Save the old contents
          !---------------------------------------------------------------------
          DO j=1,MIN(ldc(2),lda(2))
              DO i=1,MIN(ldc(1),lda(1))
                  work_field(i,j)%uc => field(i,j)%uc
                  work_field(i,j)%fc => field(i,j)%fc
                  work_field(i,j)%err => field(i,j)%err
                  work_field(i,j)%bcvalue => field(i,j)%bcvalue
              ENDDO
          ENDDO
      ENDIF

      IF (ldealloc) THEN
          !---------------------------------------------------------------------
          !  Deallocate the old contents
          !---------------------------------------------------------------------
          DO j=1,ldc(2)
              DO i=1,ldc(1)
                  IF (ASSOCIATED(field(i,j)%uc)) THEN
                      DEALLOCATE(field(i,j)%uc,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_mg_alloc_field',&
     &                       'field correction FIELD%UC',__LINE__,info)
                      ENDIF
                      NULLIFY(field(i,j)%uc)
                  ENDIF
                  IF (ASSOCIATED(field(i,j)%fc)) THEN
                      DEALLOCATE(field(i,j)%fc,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_mg_alloc_field',&
     &                       'error restriction FIELD%FC',__LINE__,info)
                      ENDIF
                      NULLIFY(field(i,j)%fc)
                  ENDIF
                  IF (ASSOCIATED(field(i,j)%err)) THEN
                      DEALLOCATE(field(i,j)%err,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_mg_alloc_field',&
     &                       'residual FIELD%ERR',__LINE__,info)
                      ENDIF
                      NULLIFY(field(i,j)%err)
                  ENDIF
                  IF (ASSOCIATED(field(i,j)%bcvalue)) THEN
                      DEALLOCATE(field(i,j)%bcvalue,STAT=info)
                      IF (info .NE. 0) THEN
                          info = ppm_error_error
                          CALL ppm_error(ppm_err_dealloc,'ppm_mg_alloc_field',&
     &                       'boundary values FIELD%BCVALUE',__LINE__,info)
                      ENDIF
                      NULLIFY(field(i,j)%bcvalue)
                  ENDIF
              ENDDO
          ENDDO
      ENDIF

      IF (lrealloc) THEN
          !---------------------------------------------------------------------
          !  Deallocate old pointer array
          !---------------------------------------------------------------------
          DEALLOCATE(field,STAT=info)
          IF (info .NE. 0) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_dealloc,'ppm_mg_alloc_field',  &
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
      CALL substop('ppm_mg_alloc_field',t0,info)
      RETURN
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_field_2d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_field_2d_sca_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_field_3d_sca_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_field_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_field_2d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_field_2d_vec_d
#endif
#elif  __MESH_DIM == __3D
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_field_3d_vec_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_alloc_field_3d_vec_d
#endif
#endif
#endif


