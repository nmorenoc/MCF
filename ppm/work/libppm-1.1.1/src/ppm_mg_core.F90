#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-----------------------------------------------------------------------
      !  Subroutine   :            ppm_mg_core    
      !-----------------------------------------------------------------------
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_mg_core.f,v $
      !  Revision 1.4  2006/05/15 14:48:57  kotsalie
      !  cosmetics
      !
      !  Revision 1.3  2004/10/14 08:45:12  kotsalie
      !
      !  Changed the criterium for going down to coarser levels
      !
      !  Revision 1.2  2004/09/23 12:40:23  kotsalie
      !  Added logical variable for distinguishing between v and w cycle
      !
      !  Revision 1.1  2004/09/22 18:31:37  kotsalie
      !  MG new version
      !
      !-----------------------------------------------------------------------  
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !------------------------------------------------------------------------ 
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      RECURSIVE SUBROUTINE ppm_mg_core_2d_sca_s(mlev,iter1,iter2,info)
#elif __KIND == __DOUBLE_PRECISION
      RECURSIVE SUBROUTINE ppm_mg_core_2d_sca_d(mlev,iter1,iter2,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      RECURSIVE SUBROUTINE ppm_mg_core_3d_sca_s(mlev,iter1,iter2,info)
#elif __KIND == __DOUBLE_PRECISION
      RECURSIVE SUBROUTINE ppm_mg_core_3d_sca_d(mlev,iter1,iter2,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      RECURSIVE SUBROUTINE ppm_mg_core_2d_vec_s(mlev,iter1,iter2,info)
#elif __KIND == __DOUBLE_PRECISION
      RECURSIVE SUBROUTINE ppm_mg_core_2d_vec_d(mlev,iter1,iter2,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      RECURSIVE SUBROUTINE ppm_mg_core_3d_vec_s(mlev,iter1,iter2,info)
#elif __KIND == __DOUBLE_PRECISION
      RECURSIVE SUBROUTINE ppm_mg_core_3d_vec_d(mlev,iter1,iter2,info)
#endif
#endif
#endif



        !-----------------------------------------------------------------------
        ! Includes
        !-----------------------------------------------------------------------
#include "ppm_define.h"

        !-----------------------------------------------------------------------
        !  Modules 
        !-----------------------------------------------------------------------
        USE ppm_module_data
        USE ppm_module_write
        USE ppm_module_data_mg
        USE ppm_module_substart 
        USE ppm_module_substop 
        USE ppm_module_mg_prolong 
        USE ppm_module_mg_restrict 
        USE ppm_module_error 
        USE ppm_module_mg_smooth
        USE ppm_module_mg_res
       

        IMPLICIT NONE

#if    __KIND == __SINGLE_PRECISION
        INTEGER, PARAMETER :: MK = ppm_kind_single
#else
        INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
        !---------------------------------------------------------------------- 
        !  Arguments     
        !----------------------------------------------------------------------

        INTEGER,                   INTENT(IN   )   ::  mlev
        INTEGER,                   INTENT(IN   )   ::  iter1
        INTEGER,                   INTENT(IN   )   ::  iter2
        INTEGER,                   INTENT(INOUT)   ::  info

        !----------------------------------------------------------------------
        !Local variables
        !----------------------------------------------------------------------
        REAL(MK)                             :: t0
        REAL(MK)                             :: scale1,scale2 
        REAL(MK)                             :: E,res
        INTEGER                              :: isub,i,j
        INTEGER                              :: ilda
        CHARACTER(LEN=256)                   :: cbuf 
        REAL(MK)                             :: c1,c2,c3,c4 
        INTEGER                              :: ncalls=0
        REAL(MK)                             :: rdx2,rdy2 
        REAL(MK)                             :: dxl,dyl 
        REAL(MK)                             :: dx,dy 
#if __MESH_DIM == __3D
        REAL(MK)                             :: c5,dzl,scale3
        REAL(MK)                             :: dz,rdz2
        INTEGER                              :: k
#endif

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



        !-----------------------------------------------------------------------
        !Externals
        !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !Initialize
        !-----------------------------------------------------------------------

        CALL substart('ppm_mg_core',t0,info)


        !---------------------------------------------------------------------  
        !  Check arguments
        !----------------------------------------------------------------------
        IF (ppm_debug .GT. 0) THEN
            IF (mlev.LE.1) THEN
                  info = ppm_error_error
                  CALL ppm_error(ppm_err_argument,'ppm_mg_core',  &
     &                'level must be >1',__LINE__,info)
                  GOTO 9999
            ENDIF
        ENDIF

        !-----------------------------------------------------------------------
        !Definition of necessary variables and allocation of arrays
        !-----------------------------------------------------------------------

#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
#if __DIM == __SFIELD
        mgfield=>mgfield_2d_sca_s
#elif __DIM == __VFIELD
        mgfield=>mgfield_2d_vec_s
#endif
        rdx2=rdx2_s
        rdy2=rdy2_s
        dx = dx_s
        dy = dy_s
#elif __KIND == __DOUBLE_PRECISION
#if __DIM == __SFIELD
        mgfield=>mgfield_2d_sca_d
#elif __DIM == __VFIELD
        mgfield=>mgfield_2d_vec_d
#endif
        rdx2=rdx2_d
        rdy2=rdy2_d
        dx = dx_d
        dy = dy_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
#if __DIM == __SFIELD
        mgfield=>mgfield_3d_sca_s
#elif __DIM == __VFIELD
        mgfield=>mgfield_3d_vec_s
#endif
        rdx2=rdx2_s
        rdy2=rdy2_s
        rdz2=rdz2_s
        dx = dx_s
        dy = dy_s
        dz = dz_s
#elif __KIND == __DOUBLE_PRECISION
#if __DIM == __SFIELD
        mgfield=>mgfield_3d_sca_d
#elif __DIM == __VFIELD
        mgfield=>mgfield_3d_vec_d
#endif
        rdx2=rdx2_d
        rdy2=rdy2_d
        rdz2=rdz2_d
        dx = dx_d
        dy = dy_d
        dz = dz_d
#endif
#endif

        !-------------------------------------------------------------------
        ! restrict the solution from the previous fine grid to the current
        ! coarser grid
        !------------------------------------------------------------------
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_restrict_2d_sca_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_restrict_2d_sca_d(mlev,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_restrict_3d_sca_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_restrict_3d_sca_d(mlev,info)
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_restrict_2d_vec_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_restrict_2d_vec_d(mlev,info)
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_restrict_3d_vec_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_restrict_3d_vec_d(mlev,info)
#endif
#endif
#endif


        !------------------------------------------------------------------
        !Initiation
        !------------------------------------------------------------------
#if __DIM == __SFIELD
#if __MESH_DIM == __2D
        scale1=REAL((factor(1)*factor(1))**(mlev-1),MK) 
        scale2=REAL((factor(2)*factor(2))**(mlev-1),MK)
        c2  =  rdx2/scale1
        c3  =  rdy2/scale2  
        c1  =  1.0_MK/(2.0_MK*(c2+c3)) 
        c4  =  1.0_MK/c1
        dxl =  dx*SQRT(scale1)
        dyl =  dy*SQRT(scale2)
          
        !--------------------------------------------------------------------
        !Compute correction using gauss-seidel algorithm
        !--------------------------------------------------------------------
        CALL ppm_mg_smooth_sca(iter1,mlev,c1,c2,c3,info) 
        !-------------------------------------------------------------------
        !Compute residual
        !-------------------------------------------------------------------
        CALL ppm_mg_res_sca(mlev,c1,c2,c3,c4,E,info) 

        !--------------------------------------------------------------------
        !Go to the next (coarser) multigrid level if the solution is
        !not converged
        !--------------------------------------------------------------------
       IF (l_print) THEN
        WRITE(cbuf,*) 'E:',E
        CALL PPM_WRITE(ppm_rank,'mg_core',cbuf,info)
       ENDIF 


        IF (mlev.LT.maxlev) THEN 
#if __KIND == __SINGLE_PRECISION

           CALL ppm_mg_core_2d_sca_s(mlev+1,iter1,iter2,info)   

        IF (w_cycle) THEN
          CALL ppm_mg_prolong_2d_sca_s(mlev,info)
          CALL ppm_mg_smooth_sca(iter2,mlev,c1,c2,c3,info)
          CALL ppm_mg_res_sca(mlev,c1,c2,c3,c4,E,info) 
          CALL ppm_mg_core_2d_sca_s(mlev+1,iter1,iter2,info)   
         ENDIF

#elif __KIND == __DOUBLE_PRECISION
           CALL ppm_mg_core_2d_sca_d(mlev+1,iter1,iter2,info)   

         IF (w_cycle) THEN
           CALL ppm_mg_prolong_2d_sca_d(mlev,info)
           CALL ppm_mg_smooth_sca(iter2,mlev,c1,c2,c3,info)
           CALL ppm_mg_res_sca(mlev,c1,c2,c3,c4,E,info) 
           CALL ppm_mg_core_2d_sca_d(mlev+1,iter1,iter2,info)   
        ENDIF

#endif
        ELSE
           GOTO 9999
        ENDIF

        !---------------------------------------------------------------------
        !else GO BACK TO A FINER LEVEL AND CONTINUE RECUSRSIVELY TO 
        !THE NEXT FINER LEVELS
        !---------------------------------------------------------------------
 

#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_prolong_2d_sca_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_prolong_2d_sca_d(mlev,info)
#endif


        scale1=REAL((factor(1)*factor(1))**(mlev-1),MK)
        scale2=REAL((factor(2)*factor(2))**(mlev-1),MK)
        c2  =  rdx2/scale1
        c3  =  rdy2/scale2
        c1  =  1.0_MK/(2.0_MK*(c2+c3))
        dxl =  dx*SQRT(scale1)
        dyl =  dy*SQRT(scale2)


        !--------------------------------------------------------------------
        !Solve for the prolongated corrections
        !--------------------------------------------------------------------

        CALL ppm_mg_smooth_sca(iter2,mlev,c1,c2,c3,info)

        !--------------------------------------------------------------------
        !Return
        !--------------------------------------------------------------------

#elif __MESH_DIM == __3D

        scale1=REAL((factor(1)*factor(1))**(mlev-1),MK) 
        scale2=REAL((factor(2)*factor(2))**(mlev-1),MK)
        scale3=REAL((factor(3)*factor(3))**(mlev-1),MK)
        c2  =  rdx2/scale1
        c3  =  rdy2/scale2  
        c4  =  rdz2/scale3  
        c1  =  1.0_MK/(2.0_MK*(c2+c3+c4)) 
        c5  =  1.0_MK/c1
        dxl =  dx*SQRT(scale1)
        dyl =  dy*SQRT(scale2)
        dzl =  dz*SQRT(scale3)


        !--------------------------------------------------------------------
        !Compute correction using gauss-seidel algorithm
       !--------------------------------------------------------------------
   
        CALL ppm_mg_smooth_sca(iter1,mlev,c1,c2,c3,c4,info) 

        !-------------------------------------------------------------------
        !Compute residual
        !-------------------------------------------------------------------
        CALL ppm_mg_res_sca(mlev,c1,c2,c3,c4,c5,E,info) 
        !--------------------------------------------------------------------
        !Go to the next (coarser) multigrid level if the solution is
        !not converged
        !--------------------------------------------------------------------
       IF (l_print) THEN
        WRITE(cbuf,*) 'E:',E
        CALL PPM_WRITE(ppm_rank,'mg_core',cbuf,info)
       ENDIF 

        IF (mlev.LT.maxlev) THEN 
           
#if __KIND == __SINGLE_PRECISION
           CALL ppm_mg_core_3d_sca_s(mlev+1,iter1,iter2,info)   
         IF (w_cycle) THEN
           CALL ppm_mg_prolong_3d_sca_s(mlev,info)
           CALL ppm_mg_smooth_sca(iter2,mlev,c1,c2,c3,c4,info)
           CALL ppm_mg_res_sca(mlev,c1,c2,c3,c4,c5,E,info) 
           CALL ppm_mg_core_3d_sca_s(mlev+1,iter1,iter2,info) 
         ENDIF  
#elif __KIND == __DOUBLE_PRECISION
           CALL ppm_mg_core_3d_sca_d(mlev+1,iter1,iter2,info)   
         IF (w_cycle) THEN
           CALL ppm_mg_prolong_3d_sca_d(mlev,info)
           CALL ppm_mg_smooth_sca(iter2,mlev,c1,c2,c3,c4,info)
           CALL ppm_mg_res_sca(mlev,c1,c2,c3,c4,c5,E,info) 
           CALL ppm_mg_core_3d_sca_d(mlev+1,iter1,iter2,info)   
         ENDIF  
#endif
        ELSE
           GOTO 9999
        ENDIF

        !---------------------------------------------------------------------
        !ELSE GO BACK TO A FINER LEVEL AND CONTINUE RECUSRSIVELY TO 
        !THE NEXT FINER LEVELS
        !---------------------------------------------------------------------

#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_prolong_3d_sca_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_prolong_3d_sca_d(mlev,info)
#endif

        scale1=REAL((factor(1)*factor(1))**(mlev-1),MK)
        scale2=REAL((factor(2)*factor(2))**(mlev-1),MK)
        scale3=REAL((factor(3)*factor(3))**(mlev-1),MK)
        c2  =  rdx2/scale1
        c3  =  rdy2/scale2
        c4  =  rdz2/scale3
        c1  =  1.0_MK/(2.0_MK*(c2+c3+c4))
        dxl =  dx*SQRT(scale1)
        dyl =  dy*SQRT(scale2)
        dzl =  dz*SQRT(scale3)


        !--------------------------------------------------------------------
        !Solve for the prolongated corrections
        !--------------------------------------------------------------------

        CALL ppm_mg_smooth_sca(iter2,mlev,c1,c2,c3,c4,info)

        !--------------------------------------------------------------------
        !--------------------------------------------------------------------

#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
        scale1=REAL((factor(1)*factor(1))**(mlev-1),MK) 
        scale2=REAL((factor(2)*factor(2))**(mlev-1),MK)
        c2  =  rdx2/scale1
        c3  =  rdy2/scale2  
        c1  =  1.0_MK/(2.0_MK*(c2+c3)) 
        c4  =  1.0_MK/c1
        dxl =  dx*SQRT(scale1)
        dyl =  dy*SQRT(scale2)
          
        !--------------------------------------------------------------------
        !Compute correction using gauss-seidel algorithm
        !--------------------------------------------------------------------
        CALL ppm_mg_smooth_vec(iter1,mlev,c1,c2,c3,info) 
        !-------------------------------------------------------------------
        !Compute residual
        !-------------------------------------------------------------------

        CALL ppm_mg_res_vec(mlev,c1,c2,c3,c4,E,info) 
        !--------------------------------------------------------------------
        !Go to the next (coarser) multigrid level if the solution is
        !not converged
        !--------------------------------------------------------------------
       IF (l_print) THEN
        WRITE(cbuf,*) 'E:',E
        CALL PPM_WRITE(ppm_rank,'mg_core',cbuf,info)
       ENDIF 

        IF (mlev.LT.maxlev) THEN 
#if __KIND == __SINGLE_PRECISION
           CALL ppm_mg_core_2d_vec_s(mlev+1,iter1,iter2,info)   
          IF (w_cycle) THEN
           CALL ppm_mg_prolong_2d_vec_s(mlev,info)
           CALL ppm_mg_smooth_vec(iter2,mlev,c1,c2,c3,info)
           CALL ppm_mg_res_vec(mlev,c1,c2,c3,c4,E,info) 
           CALL ppm_mg_core_2d_vec_s(mlev+1,iter1,iter2,info)   
          ENDIF 
#elif __KIND == __DOUBLE_PRECISION
           CALL ppm_mg_core_2d_vec_d(mlev+1,iter1,iter2,info)  
          IF (w_cycle) THEN  
           CALL ppm_mg_prolong_2d_vec_d(mlev,info)
           CALL ppm_mg_smooth_vec(iter2,mlev,c1,c2,c3,info)
           CALL ppm_mg_res_vec(mlev,c1,c2,c3,c4,E,info)
           CALL ppm_mg_core_2d_vec_d(mlev+1,iter1,iter2,info)  
          ENDIF 
#endif
        ELSE
           GOTO 9999
        ENDIF

        !---------------------------------------------------------------------
        !else GO BACK TO A FINER LEVEL AND CONTINUE RECUSRSIVELY TO 
        !THE NEXT FINER LEVELS
        !---------------------------------------------------------------------
#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_prolong_2d_vec_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_prolong_2d_vec_d(mlev,info)
#endif
        scale1=REAL((factor(1)*factor(1))**(mlev-1),MK)
        scale2=REAL((factor(2)*factor(2))**(mlev-1),MK)
        c2  =  rdx2/scale1
        c3  =  rdy2/scale2
        c1  =  1.0_MK/(2.0_MK*(c2+c3))
        dxl =  dx*SQRT(scale1)
        dyl =  dy*SQRT(scale2)


        !--------------------------------------------------------------------
        !Solve for the prolongated corrections
        !--------------------------------------------------------------------

        CALL ppm_mg_smooth_vec(iter2,mlev,c1,c2,c3,info)

        !--------------------------------------------------------------------
        !Return
        !--------------------------------------------------------------------

#elif __MESH_DIM == __3D

        scale1=REAL((factor(1)*factor(1))**(mlev-1),MK) 
        scale2=REAL((factor(2)*factor(2))**(mlev-1),MK)
        scale3=REAL((factor(3)*factor(3))**(mlev-1),MK)
        c2  =  rdx2/scale1
        c3  =  rdy2/scale2  
        c4  =  rdz2/scale3  
        c1  =  1.0_MK/(2.0_MK*(c2+c3+c4)) 
        c5  =  1.0_MK/c1
        dxl =  dx*SQRT(scale1)
        dyl =  dy*SQRT(scale2)
        dzl =  dz*SQRT(scale3)


        !--------------------------------------------------------------------
        !Compute correction using gauss-seidel algorithm
       !--------------------------------------------------------------------
   
        CALL ppm_mg_smooth_vec(iter1,mlev,c1,c2,c3,c4,info) 

        !-------------------------------------------------------------------
        !Compute residual
        !-------------------------------------------------------------------
        CALL ppm_mg_res_vec(mlev,c1,c2,c3,c4,c5,E,info) 

        !--------------------------------------------------------------------
        !Go to the next (coarser) multigrid level if the solution is
        !not converged
        !--------------------------------------------------------------------
       IF (l_print) THEN
        WRITE(cbuf,*) 'E:',E
        CALL PPM_WRITE(ppm_rank,'mg_core',cbuf,info)
       ENDIF 
        IF (mlev.LT.maxlev) THEN 
           
#if __KIND == __SINGLE_PRECISION
           CALL ppm_mg_core_3d_vec_s(mlev+1,iter1,iter2,info)   
          IF (w_cycle) THEN
           CALL ppm_mg_prolong_3d_vec_s(mlev,info)
           CALL ppm_mg_smooth_vec(iter2,mlev,c1,c2,c3,c4,info)
           CALL ppm_mg_res_vec(mlev,c1,c2,c3,c4,c5,E,info) 
           CALL ppm_mg_core_3d_vec_s(mlev+1,iter1,iter2,info)
          ENDIF   
#elif __KIND == __DOUBLE_PRECISION
           CALL ppm_mg_core_3d_vec_d(mlev+1,iter1,iter2,info)   
           IF (w_cycle) THEN
           CALL ppm_mg_prolong_3d_vec_d(mlev,info)
           CALL ppm_mg_smooth_vec(iter2,mlev,c1,c2,c3,c4,info)
           CALL ppm_mg_res_vec(mlev,c1,c2,c3,c4,c5,E,info) 
           CALL ppm_mg_core_3d_vec_d(mlev+1,iter1,iter2,info) 
           ENDIF  
#endif
        ELSE
           GOTO 9999
        ENDIF

        !---------------------------------------------------------------------
        !ELSE GO BACK TO A FINER LEVEL AND CONTINUE RECUSRSIVELY TO 
        !THE NEXT FINER LEVELS
        !---------------------------------------------------------------------

#if __KIND == __SINGLE_PRECISION
        CALL ppm_mg_prolong_3d_vec_s(mlev,info)
#elif __KIND == __DOUBLE_PRECISION
        CALL ppm_mg_prolong_3d_vec_d(mlev,info)
#endif

        scale1=REAL((factor(1)*factor(1))**(mlev-1),MK)
        scale2=REAL((factor(2)*factor(2))**(mlev-1),MK)
        scale3=REAL((factor(3)*factor(3))**(mlev-1),MK)
        c2  =  rdx2/scale1
        c3  =  rdy2/scale2
        c4  =  rdz2/scale3
        c1  =  1.0_MK/(2.0_MK*(c2+c3+c4))
        dxl =  dx*SQRT(scale1)
        dyl =  dy*SQRT(scale2)
        dzl =  dz*SQRT(scale3)


        !--------------------------------------------------------------------
        !Solve for the prolongated corrections
        !--------------------------------------------------------------------

        CALL ppm_mg_smooth_vec(iter2,mlev,c1,c2,c3,c4,info)

        !--------------------------------------------------------------------
        !--------------------------------------------------------------------

#endif
#endif


9999    CONTINUE
        CALL substop('ppm_mg_core',t0,info)
        RETURN

#if __DIM == __SFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_core_2d_sca_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_core_2d_sca_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_core_3d_sca_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_core_3d_sca_d
#endif
#endif
#elif __DIM == __VFIELD
#if __MESH_DIM == __2D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_core_2d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_core_2d_vec_d
#endif
#elif __MESH_DIM == __3D
#if __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_mg_core_3d_vec_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_mg_core_3d_vec_d
#endif
#endif
#endif


