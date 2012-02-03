program TestRunner
  use ppm_module_init_test
  use ppm_module_mktopo_test
  use ppm_module_map_part_ghost_test
  use ppm_module_get_revision_test
  implicit none

  integer          :: i
  character(len=3) :: arg

  print *, 'staring tests runner'

  if (iargc() < 1) then 
     print *, __FILE__, __LINE__, "TestRunner should have an argumetn"
     stop 20
  endif

  call getarg(1, arg)
  read (arg, '(i3)') i

  select case(i) 
  case(1)
     call ppm_module_init_run
  case(2)
     call ppm_module_mktopo_run
  case(3)
     call ppm_module_map_part_ghost_run
  case(4)
     call ppm_module_get_revision_run
  case default
     print *, __FILE__, __LINE__, "Unknown test case"
  end select

  print *, 'end tests'
end program TestRunner
