#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_get_revision
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Return version and revision numbers of ppm as well
      !                 as revision date from RCS/CVS.
      !
      !  Input        :
      !
      !  Input/output : 
      !
      !  Output       : version       (C) ppm library version (release)
      !                 revision      (C) CVS revision and revision data
      !                 info          (I) error status. 0 on success.
      !
      !  Remarks      : All information is taken from the CVS/RCS system.
      !                 No hard-coded values need to be updated in this
      !                 file. 
      !
      !                 For the information to be up to date, commit this
      !                 file whenever a new release is done.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_get_revision.f,v $
      !  Revision 1.8  2004/11/12 15:30:00  ivos
      !  New revision.
      !
      !  Revision 1.7  2004/10/01 16:33:33  ivos
      !  cosmetics.
      !
      !  Revision 1.6  2004/10/01 16:08:59  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.5  2004/07/26 07:45:29  ivos
      !  Updated to use single-interface modules. Adapted all USE statements.
      !
      !  Revision 1.4  2004/07/16 14:46:25  ivos
      !  Added check for ppm_initialized.
      !
      !  Revision 1.3  2004/06/01 13:24:57  ivos
      !  bugfix: forgot to USE ppm_module_util. Funny, pgf90 did not complain...
      !
      !  Revision 1.2  2004/06/01 11:31:18  ivos
      !  bugfix: added initialization for idx1 and idx2.
      !
      !  Revision 1.1  2004/06/01 09:25:43  ivos
      !  Initial implementation.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_get_revision(version,revision,info)
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_substart
      USE ppm_module_substop
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(  OUT) :: version,revision
      INTEGER         , INTENT(  OUT) :: info

      version = PACKAGE_VERSION
      revision = 'undefined'
      ! never fails
      info = 0

      RETURN
      END SUBROUTINE ppm_get_revision
