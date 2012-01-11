AM_FCFLAGS   = -O2 -D__MPI -D__FLOW_DEVELOPED -D__COLLOID_SEPARATE_FILE
bin_PROGRAMS = mcf
mcf_SOURCES  = mcf.F90\
	       mcf_header.F90\
	       Class_Tool.F90\
	       Class_Random.F90\
	       Class_Debug.F90\
	       Class_Control.F90\
               Class_Boundary.F90\
               Class_Colloid.F90\
	       Class_Physics.F90\
               Class_Rhs.F90\
	       Class_IO.F90\
	       Class_StateEquation.F90\
	       Class_Kernel.F90\
               Class_Technique.F90\
	       Class_Particles.F90\
	       Class_Statistic.F90\
	       Class_Marching.F90

Makefile.dep:
	$(MAKEDEPF90) $(mcf_SOURCES) > $@

clean-local:
	-rm *.mod Makefile.dep

-include Makefile.dep



