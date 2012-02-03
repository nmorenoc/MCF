      !-------------------------------------------------------------------------
      !  Module       :                    ppm_define
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Describes and sets the compile-time DEFINEs
      !
      !  Remarks      : USE_MPI     Set to build the parallel version. Serial i
      !                           version is built if not set.
      !                 __VECTOR  Set to build version which is optimized for 
      !                           vector processors.
      !                 __SXF90   Enable NEC SX Fortran vector directives
      !                 __XLF     Enable IBM Fortran directives
      !                 __POWERPC Enable PowerPC-specific optimizations
      !                 __NOMICROINSTRUCTIONS Do interpolations without
      !                           f90 microinstructions (colon notation). 
      !                 __ETIME   Set if the operating system provides the 
      !                           etime timing facility. If not set, CPU_TIME 
      !                           is used.
      !                 __Linux   Set if OS is some kind of Linux.
      !                 HAVE_LIBFFTW3    Set to embed FFTW-library
      !                 __MATHKEISAN   Set to embed MathKeisan Library
      !                 HAVE_LIBMETIS   Set to compile ppm with support for the
      !                           METIS library for mesh decomposition and 
      !                           assignment of subdomains to processors
      !                 __CRAYFISHPACK  Set to compile ppm with the 
      !                                 crayfishpack field solver
      !                 __HYPRE   Set to compile ppm with support for 
      !                           algebraic multigrid solvers from Hypre.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_define.h,v $
      !  Revision 1.14  2006/09/04 18:34:42  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.13  2005/06/23 18:24:27  ivos
      !  Added __POWERPC and __NOMICROINSTRUCTIONS
      !
      !  Revision 1.12  2004/11/05 18:16:30  ivos
      !  Updated header comment.
      !
      !  Revision 1.11  2004/11/05 18:13:16  michaebe
      !  added __XLF
      !
      !  Revision 1.10  2004/11/03 11:11:26  hiebers
      !  introduced new variable __MATHKEISAN
      !
      !  Revision 1.9  2004/06/15 08:40:58  ivos
      !  Added SXF90 define.
      !
      !  Revision 1.8  2004/02/25 14:03:57  walther
      !  Let us try to have the commented defines as in this version.
      !
      !  Revision 1.7  2004/02/24 10:36:11  hiebers
      !  merged
      !
      !  Revision 1.6  2004/02/12 17:08:21  ivos
      !  update header comments.
      !
      !  Revision 1.5  2004/02/11 16:23:07  ivos
      !  Added METIS, HYPRE and CRAYFISHPACK defines.
      !
      !  Revision 1.4  2004/02/10 17:10:47  hiebers
      !  add comment
      !
      !  Revision 1.3  2004/02/10 16:32:14  hiebers
      !  added define HAVE_LIBFFTW3
      !
      !  Revision 1.2  2004/01/23 17:25:16  ivos
      !  Update.
      !
      !  Revision 1.1  2004/01/13 12:18:43  ivos
      !  New implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Define 
      !-------------------------------------------------------------------------
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define __Linux

/* 
#define __Linux
#define __ETIME
#define __VECTOR 
#define __MATHKEISAN
#define __SXF90
#define __CRAYFISHPACK 
#define __HYPRE 
#define __NOMICROINSTRUCTIONS
*/
