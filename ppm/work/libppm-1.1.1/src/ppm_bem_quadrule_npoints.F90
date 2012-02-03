#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_bem_quadrule_npoints
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Returns the number of quadrature points for the
      !                 given rule
      !
      !  Input        : rule    (I) : the quadrature rule
      !
      !  Input/output : 
      !
      !  Output       : nqp     (I) : number of quadrature points
      !                 info    (I) : return status. 0 on success.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_bem_quadrule_npoints.f,v $
      !  Revision 1.8  2006/09/04 18:34:41  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.7  2006/05/11 10:44:07  pchatela
      !  Added Hammer rules
      !  Switch loop rather than if then else
      !
      !  Revision 1.6  2004/10/01 16:33:31  ivos
      !  cosmetics.
      !
      !  Revision 1.5  2004/10/01 16:08:56  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.4  2004/07/27 07:06:48  ivos
      !  Renamed parameters after moving their definition to ppm_param.h.
      !
      !  Revision 1.3  2004/07/26 11:46:53  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.2  2004/07/26 07:46:36  ivos
      !  Changed to use single-interface modules. Updated all USE statements.
      !
      !  Revision 1.1  2004/07/16 08:33:39  oingo
      !  Initial release. Not tested
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_bem_quadrule_npoints(rule,nqp,info)

      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN   ) :: rule
      INTEGER, INTENT(  OUT) :: nqp               
      INTEGER, INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)  :: t0
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_bem_quadrule_npoints',t0,info)

      !-------------------------------------------------------------------------
      !  Check for the rule and set the number of points accordingly
      !-------------------------------------------------------------------------
      SELECT CASE (rule)
      CASE(ppm_param_bem_quadrule_center)
         nqp = 1
      CASE(ppm_param_bem_quadrule_nodes)
         nqp = 3
      CASE(ppm_param_bem_quadrule_edges)
         nqp = 3
      CASE(ppm_param_bem_quadrule_cne)
         nqp = 7
      CASE(ppm_param_bem_quadrule_stroud,ppm_param_bem_quadrule_hammer7)
         nqp = 7
      CASE(ppm_param_bem_quadrule_hammer3)
         nqp = 3
      CASE(ppm_param_bem_quadrule_hammer4)
         nqp = 4
      CASE(ppm_param_bem_quadrule_hammer12)
         nqp = 12
      CASE DEFAULT
         nqp = 0
      END SELECT
      
      IF (nqp .EQ. 0) THEN
         info = ppm_error_error
         CALL ppm_error(ppm_err_argument,'ppm_bem_quadrule_npoints',   &
     &        'Unknown quadrature rule',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_bem_quadrule_npoints',t0,info)
      RETURN

      END SUBROUTINE ppm_bem_quadrule_npoints

