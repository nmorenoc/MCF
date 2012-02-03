#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module ppm_module_get_revision_test
  implicit none
contains
  subroutine ppm_module_get_revision_run
    use ppm_module_get_revision, only : ppm_get_revision
    use mod_unit
    implicit none
    character(len=100) :: version, revision
    integer            :: info

    call unit_init(1)
    call ppm_get_revision(version,revision,info)
    call unit_assert_equal('test version', version, PACKAGE_VERSION)
    call unit_finalize()

  end subroutine ppm_module_get_revision_run
end module ppm_module_get_revision_test
