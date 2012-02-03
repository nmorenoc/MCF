#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                ppm_fdsolver_map_2d
      !-------------------------------------------------------------------------
      !
      !  Purpose      : maps the data from topology topo_ids(1)) / mesh 
      !               (mesh_ids(1)) to topology (topo_ids(2))/mesh (mesh_ids(2))
      !                 
      !
      !  Input        :  
      !                 topo_ids(2)        (I)     
      !                                         first: current topology   
      !                                         second: destination topology   
      !                 mesh_ids(2)        (I)
      !                                         first: current mesh of the data
      !                                         second: destination mesh   
      !                                
      !
      !  Input/output : data_fv(:,:,:,:)   (F)  data to be mapped
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
      !  $Log: ppm_fdsolver_map_2d.f,v $
      !  Revision 1.9  2006/09/04 18:34:43  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.8  2005/05/03 13:33:41  hiebers
      !  Bugfix: call ppm_check_meshid with internal topo ID
      !
      !  Revision 1.7  2004/10/01 16:08:57  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.6  2004/08/31 15:14:03  hiebers
      !  added argument check for topo and mesh ids
      !
      !  Revision 1.5  2004/08/19 13:14:54  hiebers
      !  debugged scalar/vector version
      !
      !  Revision 1.3  2004/07/26 13:49:16  ivos
      !  Removed Routines sections from the header comment.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_map_2d_sca_s(data_fv, topo_ids, mesh_ids, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_map_2d_sca_d(data_fv, topo_ids, mesh_ids, info)
#elif __KIND == __COMPLEX
      SUBROUTINE ppm_fdsolver_map_2d_sca_c(data_fv, topo_ids, mesh_ids, info)
#elif __KIND == __DOUBLE_COMPLEX
      SUBROUTINE ppm_fdsolver_map_2d_sca_cc(data_fv, topo_ids, mesh_ids,info)
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_fdsolver_map_2d_vec_s(data_fv,lda,topo_ids, mesh_ids, info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_fdsolver_map_2d_vec_d(data_fv,lda,topo_ids, mesh_ids, info)
#elif __KIND == __COMPLEX
      SUBROUTINE ppm_fdsolver_map_2d_vec_c(data_fv,lda,topo_ids, mesh_ids, info)
#elif __KIND == __DOUBLE_COMPLEX
      SUBROUTINE ppm_fdsolver_map_2d_vec_cc(data_fv,lda,topo_ids, mesh_ids,info)
#endif
#endif



      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_write
      USE ppm_module_error
      USE ppm_module_map_field
      USE ppm_module_check_topoid
      USE ppm_module_check_meshid



      IMPLICIT NONE

#if   __KIND == __SINGLE_PRECISION | __KIND == __COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION | __KIND == __DOUBLE_COMPLEX 
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif



      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      ! POINTER to data

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:,:),          POINTER   :: data_fv
#elif __KIND == __COMPLEX | __KIND == __DOUBLE_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:),       POINTER   :: data_fv
#endif
#endif


#if __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION | __KIND == __DOUBLE_PRECISION
      REAL(MK), DIMENSION(:,:,:,:),          POINTER   :: data_fv
#elif __KIND == __COMPLEX | __KIND == __DOUBLE_COMPLEX
      COMPLEX(MK), DIMENSION(:,:,:,:),       POINTER   :: data_fv
#endif
#endif


#if   __DIM == __VFIELD
      INTEGER,                          INTENT(IN)     :: lda
#endif  
   
      INTEGER, DIMENSION(2),            INTENT(IN)     :: topo_ids
      INTEGER, DIMENSION(2),            INTENT(IN)     :: mesh_ids
      INTEGER              ,            INTENT(  OUT)  :: info


      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! timer
      REAL(MK)                                :: t0
      ! counters
      INTEGER                                 :: k, i, j
      !Size of the data_in 
      INTEGER                                 :: from_topo, to_topo
      INTEGER                                 :: from_mesh, to_mesh
      LOGICAL                                 :: valid

#if   __DIM == __SFIELD
      INTEGER, PARAMETER                      :: lda = 1
#endif   


      INTEGER , DIMENSION(2  )                :: ghostsize
      INTEGER                                 :: maptype
      CHARACTER(LEN=ppm_char)                 :: mesg

      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_fdsolver_map_2d',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
          IF (ppm_internal_topoid(topo_ids(1)) .NE. ppm_field_topoid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_fdsolver_map_2d',  &
     &            'Passed topology is not the current topology',__LINE__,info)
              GOTO 9999
          ENDIF

          IF (topo_ids(2).GE.0) THEN
              CALL ppm_check_topoid(ppm_param_id_user,topo_ids(2),valid,info)
              IF (.NOT. valid) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_fdsolver_map_2d',  &
     &                'Topology ID (from_topo) is invalid!',__LINE__,info)
                  GOTO 9999
              ENDIF           
          ENDIF
          IF (mesh_ids(1) .GT. 0) THEN
                CALL ppm_check_meshid(ppm_param_id_user,mesh_ids(1),     &
     &              ppm_internal_topoid(topo_ids(1)),valid,info)
                IF (.NOT. valid) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_fdsolver_map_2d',  &
     &                  'Mesh ID (from_mesh) invalid!',__LINE__,info)
                    GOTO 9999
                ENDIF
          ENDIF
          IF (mesh_ids(2) .GT. 0) THEN
                CALL ppm_check_meshid(ppm_param_id_user,mesh_ids(2),     &
     &              ppm_internal_topoid(topo_ids(2)),valid,info)
                IF (.NOT. valid) THEN
                    info = ppm_error_error
                    CALL ppm_error(ppm_err_argument,'ppm_fdsolver_map_2d',  &
     &                  'Mesh ID (to_mesh) invalid!',__LINE__,info)
                    GOTO 9999
                ENDIF
          ENDIF



      ENDIF
      

      !-------------------------------------------------------------------------
      ! Define source and destination 
      !-------------------------------------------------------------------------
      from_topo   = topo_ids(1)
      to_topo   = topo_ids(2)
      from_mesh = mesh_ids(1)
      to_mesh   = mesh_ids(2)
      ghostsize(1)      = 0
      ghostsize(2)      = 0

      IF (ppm_debug .GT. 0) THEN

         WRITE(mesg,'(A,I5,A,I5)' )'Mapping from topo ',from_topo,        &
    &          ', mesh ',from_mesh 
         CALL ppm_write(ppm_rank,'ppm_fdsolver_map',mesg,j)      

         WRITE(mesg,'(A,I5,A,I5)' )'          to topo ',to_topo,        &
   &             ', mesh ',to_mesh 
         CALL ppm_write(ppm_rank,'ppm_fdsolver_map',mesg,j)      
      ENDIF
 

      !-------------------------------------------------------------------------
      !  Map fields
      !-------------------------------------------------------------------------

#if   __DIM == __SFIELD
      maptype = ppm_param_map_global
      CALL ppm_map_field(DATA_fv,to_topo,from_mesh,to_mesh,ghostsize,maptype,  &
     &           info)
      maptype = ppm_param_map_push
      CALL ppm_map_field(DATA_fv,to_topo,from_mesh,to_mesh,ghostsize,maptype,  &
     &           info)
      maptype = ppm_param_map_send
      CALL ppm_map_field(DATA_fv,to_topo,from_mesh,to_mesh,ghostsize,maptype,  &
     &           info)
      maptype = ppm_param_map_pop
      CALL ppm_map_field(DATA_fv,to_topo,from_mesh,to_mesh,ghostsize,maptype,  &
     &           info)

#elif __DIM == __VFIELD     
      maptype = ppm_param_map_global;
      CALL ppm_map_field(DATA_fv,lda,to_topo,from_mesh,to_mesh,ghostsize,      &
     &           maptype,info)
      maptype = ppm_param_map_push;
      CALL ppm_map_field(DATA_fv,lda,to_topo,from_mesh,to_mesh,ghostsize,      &
     &           maptype,info)
      maptype = ppm_param_map_send;
      CALL ppm_map_field(DATA_fv,lda,to_topo,from_mesh,to_mesh,ghostsize,      & 
     &           maptype,info)
      maptype = ppm_param_map_pop;
      CALL ppm_map_field(DATA_fv,lda,to_topo,from_mesh,to_mesh,ghostsize,      &
     &           maptype,info)
#endif



      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_fdsolver_map_2d',t0,info)

      RETURN

#if   __DIM == __SFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_map_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_map_2d_sca_d
#elif __KIND == __COMPLEX
      END SUBROUTINE ppm_fdsolver_map_2d_sca_c
#elif __KIND == __DOUBLE_COMPLEX
      END SUBROUTINE ppm_fdsolver_map_2d_sca_cc
#endif
#endif

#if   __DIM == __VFIELD
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_fdsolver_map_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_fdsolver_map_2d_vec_d
#elif __KIND == __COMPLEX
      END SUBROUTINE ppm_fdsolver_map_2d_vec_c
#elif __KIND == __DOUBLE_COMPLEX
      END SUBROUTINE ppm_fdsolver_map_2d_vec_cc
#endif
#endif
