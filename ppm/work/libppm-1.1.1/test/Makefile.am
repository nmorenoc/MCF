bin_PROGRAMS=TestRunner
TestRunner_LDADD=${top_srcdir}/src/libppm.a

if HAVE_MPI

else 
TestRunner_LDADD+=${top_srcdir}/src/libfakempi.a
endif

TestRunner_SOURCES=TestRunner.F90 mod_unit.F90 ppm_module_init_test.F90 \
	ppm_module_mktopo_test.F90 ppm_module_map_part_ghost_test.F90 \
	ppm_module_get_revision_test.F90

AM_FCFLAGS=-I${top_srcdir}/src

-include Makefile.dep 

MAKEDEPF90_FLAGS=-DHAVE_CONFIG_H
Makefile.dep:
	$(MAKEDEPF90) $(MAKEDEPF90_FLAGS) -I ${top_builddir} $(TestRunner_SOURCES) | \
	sed -e 's,config.h,${CONFIG_HEADER},g' > "$@"



clean-local:
	-rm  *.mod