#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_mg_finalize
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine deallocates all the arrays 
      !                 
      !
      !  Input        :
      !
      !  Input/output : 
      !
      !  Output       : info    (I) 0 on success.
      !
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_mg_finalize.f,v $
      !  Revision 1.5  2006/07/21 11:30:55  kotsalie
      !  FRIDAY
      !
      !  Revision 1.3  2004/10/01 16:33:39  ivos
      !  cosmetics.
      !
      !  Revision 1.2  2004/10/01 16:09:11  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.1  2004/09/23 09:56:40  kotsalie
      !  Mg new version
      !
      !  Revision 1.4  2004/07/26 13:49:19  ivos
      !  Removed Routines sections from the header comment.
      !
      !  Revision 1.3  2004/07/26 11:59:40  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.2  2004/07/26 07:42:39  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.1  2004/06/29 14:36:49  kotsalie
      !  Commiting multigrid for further use
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
     SUBROUTINE ppm_mg_finalize_2d_sca_s(info)
#elif __KIND == __DOUBLE_PRECISION
     SUBROUTINE ppm_mg_finalize_2d_sca_d(info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
     SUBROUTINE ppm_mg_finalize_3d_sca_s(info)
#elif __KIND == __DOUBLE_PRECISION
     SUBROUTINE ppm_mg_finalize_3d_sca_d(info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
     SUBROUTINE ppm_mg_finalize_2d_vec_s(info)
#elif __KIND == __DOUBLE_PRECISION
     SUBROUTINE ppm_mg_finalize_2d_vec_d(info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
     SUBROUTINE ppm_mg_finalize_3d_vec_s(info)
#elif __KIND == __DOUBLE_PRECISION
     SUBROUTINE ppm_mg_finalize_3d_vec_d(info)
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
      USE ppm_module_data_mg
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_mg_alloc
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION | __KIND == __SINGLE_PRECISION_COMPLEX
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(INOUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(1) :: lda1 
      INTEGER, DIMENSION(2) :: lda2 
      INTEGER, DIMENSION(3) :: lda3 
      INTEGER, DIMENSION(4) :: lda4 
      INTEGER, DIMENSION(5) :: lda5 
      INTEGER               :: iopt,i
      INTEGER               :: istat,j
      REAL(MK)              :: t0
      CHARACTER(LEN=ppm_char) :: mesg
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      TYPE(mg_field_2d_sca_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_2d_sca_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      TYPE(mg_field_3d_sca_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_3d_sca_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      TYPE(mg_field_2d_vec_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_2d_vec_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      TYPE(mg_field_3d_vec_s),DIMENSION(:,:),POINTER :: mgfield
#elif __KIND == __DOUBLE_PRECISION
      TYPE(mg_field_3d_vec_d),DIMENSION(:,:),POINTER :: mgfield
#endif
#endif
#endif

      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_mg_finalize',t0,info)
      lda1(1)=0
      lda2(1:2) = 0
      lda3(1:3) = 0
      lda4(1:4) = 0
      lda5(1:5) = 0

      !-------------------------------------------------------------------------
      !Definition of necessary variables and allocation of arrays
      !-------------------------------------------------------------------------
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      mgfield=>mgfield_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
      mgfield=>mgfield_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
       mgfield=>mgfield_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
       mgfield=>mgfield_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      mgfield=>mgfield_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      mgfield=>mgfield_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
       mgfield=>mgfield_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
       mgfield=>mgfield_3d_vec_d
#endif
#endif
#endif


      !-------------------------------------------------------------------------
      !  Deallocate global arrays (from ppm_module_multigrid)
      !-------------------------------------------------------------------------
      istat = 0
      iopt = ppm_param_dealloc
      
       CALL ppm_alloc(start,lda3,iopt,info)
       istat=istat+info  
       CALL ppm_alloc(istop,lda3,iopt,info)
       istat=istat+info  
       CALL ppm_alloc(lboundary,lda2,iopt,info)
       istat=istat+info  
       CALL ppm_alloc(max_node,lda2,iopt,info)
       istat=istat+info  
#if __DIM == __SFIELD
       CALL ppm_alloc(bcdef_sca,lda1,iopt,info)
       istat=istat+info  
#elif __DIM == __VFIELD
       CALL ppm_alloc(bcdef_vec,lda1,iopt,info)
       istat=istat+info
#endif
       CALL ppm_alloc(ghostsize,lda1,iopt,info)
       istat=istat+info  
       CALL ppm_alloc(factor,lda1,iopt,info)
       istat=istat+info  
       CALL ppm_alloc(mesh_id_g,lda1,iopt,info)
       istat=istat+info  
       CALL ppm_alloc(meshid_g,lda1,iopt,info)
       istat=istat+info  
       CALL ppm_mg_alloc(mgfield,lda2,iopt,info)
       istat = istat +info

      IF (istat .NE. 0) THEN
          WRITE(mesg,'(A,I3,A)') 'for ',istat,' mgr arrays.Pble memory leak.'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_mg_finalize',mesg,__LINE__,&
     &                                                                 info)
          GOTO 9999
      ENDIF



      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_mg_finalize',t0,info)
      RETURN

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_mg_finalize_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_mg_finalize_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_mg_finalize_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_mg_finalize_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_mg_finalize_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_mg_finalize_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
     END SUBROUTINE ppm_mg_finalize_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
     END SUBROUTINE ppm_mg_finalize_3d_vec_d
#endif
#endif
#endif

