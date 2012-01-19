# prefix to install libraries locally
# libarieris go to $PREFIX/lib 
# binaries go to  $PREFIX/bin
# headers go to $PREFIX/bin 
PREFIX=$(PWD)/local_mpich_ifort121_o0_g
BUILD_MAKE_FLAGS='-j 8'

# a directory to store tarbars
DISTDIR=$(PWD)/distfiles

# c++ compiler 
CXX=icc

# c  compiler
CC=icc

# fortran compiler 
FC=ifort
FCFLAGS=-O2
MPIFC=mpif90

# for g95
#FCFLAGS=-fno-second-underscore

# mpi wrapper to fortran compiler

# ppm
PPM_VERSION=1.1.1
PPM_DISTR_FILE=libppm-$(PPM_VERSION).tar.gz

# a place to fetch fftw library
FFTW_SITE=ftp://ftp.netbsd.org/pub/NetBSD/packages/distfiles
FFTW_VERSION=3.1.3

# a place to fetch mpich library
MPICH_VERSION=1.0.8
MPICH_SITE=http://www.mcs.anl.gov/research/projects/mpich2/downloads/tarballs/

MAKEDEPF90_SITE=http://personal.inet.fi/private/erikedelmann/makedepf90/
MAKEDEPF90=$(PREFIX)/bin/makedepf90

# a preprocessor for  ppm
CPP=cpp
RSH=/usr/bin/ssh


EXTRACT_COMAND=scripts/extract
WRKDIR=$(PWD)/work_mpich_ifort121_o0_g


all: 
	@echo "!!!!!!!!!!"
	@echo "There is no default target"
	@echo "set FC, PREFIX and DISTDIR in Makefile"
	@echo "Please refer README file for possible targets"

# compile and install 
install_all_but_mpi: install_fftw \
	install_vizing install_metis \
	install_makedepf90 install_ppm	


# install mpich localy
install_mpi: mpich_extract mpich_configure mpich_compile mpich_install 

install_all: install_mpi install_all_but_mpi

install_metis: metis_extract metis_configure metis_compile metis_install
install_vizing: vizing_extract vizing_configure vizing_compile vizing_install
install_ppm: ppm_extract ppm_configure ppm_compile ppm_install

install_fftw: fftw_extract fftw_compile

install_pse: ppm-pse_extract  ppm-pse_compile 

install_sph: ppm-sph_extract ppm-sph_compile 

install_makedepf90: makedepf90_extract makedepf90_configure makedepf90_compile makedepf90_install




# get all 'tarbars'
fetch_all: mpich_fetch fftw_fetch


mpich_fetch: 
	cd $(DISTDIR); \
	ftp $(MPICH_SITE)/$(MPICH_VERSION)/mpich2-$(MPICH_VERSION).tar.gz

mpich_extract: $(DISTDIR)/mpich2-$(MPICH_VERSION).tar.gz
	${EXTRACT_COMAND} -d ${WRKDIR} $<

mpich_configure:
	cd ${WRKDIR}/mpich2-$(MPICH_VERSION) ; \
	./configure  RSHCOMMAND="$(RSH)" MPICH2LIB_FCFLAGS=-O2  \
	FCFLAGS="$(FCFLAGS)" FC="$(FC)" F77="$(FC)" F90="$(FC)" \
	CC=$(CC) \
        --disable-romio --enable-f90 --disable-cxx --enable-fast=O3 \
        --prefix=$(PREFIX)

mpich_compile:
	cd ${WRKDIR}/mpich2-$(MPICH_VERSION) ; \
	make ${BUILD_MAKE_FLAGS}

mpich_install:
	cd ${WRKDIR}/mpich2-$(MPICH_VERSION) ; \
	make install


mpich_clean:
	-rm -rf ${WRKDIR}/mpich2-$(MPICH_VERSION)




fftw_fetch: 
	cd $(DISTDIR); \
	ftp $(FFTW_SITE)/fftw-$(FFTW_VERSION).tar.gz

fftw_extract: $(DISTDIR)/fftw-$(FFTW_VERSION).tar.gz
	${EXTRACT_COMAND} -d ${WRKDIR} $<

fftw_compile: 
	cd ${WRKDIR}/fftw-$(FFTW_VERSION); \
	./configure --prefix=$(PREFIX) CC=$(CC) FC=$(FC) F77=$(FC) F90=$(FC) FCFLAGS="$(FCFLAGS)" \
	--enable-fortran --enable-float; \
	make ${BUILD_MAKE_FLAGS} && make install; \
	./configure --prefix=$(PREFIX) CC=$(CC) FC=$(FC) F77=$(FC) F90=$(FC) FCFLAGS="$(FCFLAGS)" \
	--enable-fortran; \
	make ${BUILD_MAKE_FLAGS} && make install

fftw_clean:
	-rm -rf ${WRKDIR}/fftw-$(FFTW_VERSION)


ppm_extract: $(DISTDIR)/$(PPM_DISTR_FILE)
	${EXTRACT_COMAND} -d ${WRKDIR} $<

ppm_configure:
	cd ${WRKDIR}/libppm-$(PPM_VERSION); \
	./configure --prefix=$(PREFIX) --enable-mpi \
	MPIFC=$(MPIFC) LDFLAGS=-L$(PREFIX)/lib \
	MAKEDEPF90=$(MAKEDEPF90) FCFLAGS="-I$(PREFIX)/include -g -O0"

ppm_compile:
	cd ${WRKDIR}/libppm-$(PPM_VERSION); \
	make ${BUILD_MAKE_FLAGS}

ppm_install:
	cd ${WRKDIR}/libppm-$(PPM_VERSION); \
	make install

# ppm_patch:
# 	patch ${WRKDIR}/ppm-*/src/ppm_mg_smooth_coarse.f90 patches/aa; \
# 	patch ${WRKDIR}/ppm-*/src/ppm_mg_smooth_fine.f90 patches/ab; \
# 	patch ${WRKDIR}/ppm-*/src/ppm_tree_boxcut.f90 patches/ac

# ppm_clean: 
# 	-rm -rf ${WRKDIR}/ppm-*


vizing_extract: $(DISTDIR)/vizing-0.1.tar.gz
	${EXTRACT_COMAND} -d ${WRKDIR} $<

vizing_configure:
	cd ${WRKDIR}/vizing-*; \
	./configure --prefix=${PREFIX} CXX=$(CXX)

vizing_compile:
	cd ${WRKDIR}/vizing-*; \
	make ${BUILD_MAKE_FLAGS}

vizing_install:
	cd ${WRKDIR}/vizing-*; \
	make install

makedepf90_fetch: 
	cd $(DISTDIR); \
	ftp $(MAKEDEPF90_SITE)/makedepf90-2.8.8.tar.gz

makedepf90_extract: $(DISTDIR)/makedepf90-2.8.8.tar.gz
	${EXTRACT_COMAND} -d ${WRKDIR} $<

makedepf90_configure:
	cd ${WRKDIR}/makedepf90-*; \
	./configure --prefix=${PREFIX} CC=$(CC)

makedepf90_compile:
	cd ${WRKDIR}/makedepf90-*; \
	make ${BUILD_MAKE_FLAGS}

makedepf90_install:
	cd ${WRKDIR}/makedepf90-*; \
	make install


metis_extract: $(DISTDIR)/metis-4.1.1.tar.gz
	${EXTRACT_COMAND} -d ${WRKDIR} $<

metis_configure:
	cd ${WRKDIR}/metis-*; \
	./configure --prefix=${PREFIX} CC=$(CC)

metis_compile:
	cd ${WRKDIR}/metis-*; \
	make ${BUILD_MAKE_FLAGS}

metis_install:
	cd ${WRKDIR}/metis-*; \
	make install


ppm-pse_extract: $(DISTDIR)/ppm-pse.tar.gz
	${EXTRACT_COMAND} -d ${WRKDIR} $<


ppm-pse_compile:
	cd ${WRKDIR}/ppm-pse; \
	LD=$(MPIFC)  FC=$(MPIFC) PREFIX=$(PREFIX) FCFLAGS="$(FCFLAGS)" make  

ppm-sph_extract: $(DISTDIR)/ppm-sph.tar.gz
	${EXTRACT_COMAND} -d ${WRKDIR} $<


ppm-sph_compile:
	cd ${WRKDIR}/ppm-sph; \
	FC=$(MPIFC) FCFLAGS=$(FCFLAGS) LD=$(MPIFC) PREFIX=$(PREFIX) make  


ppm-dpd_extract: $(DISTDIR)/ppm-dpd.tar.gz
	${EXTRACT_COMAND} -d ${WRKDIR} $< 

ppm-dpd_compile:
	cd ${WRKDIR}/ppm-dpd; \
	LD=$(MPIFC)  FC=$(MPIFC) PREFIX=$(PREFIX) LDFLAGS=-L$(PREFIX)/lib FCFLAGS="$(FCFLAGS)" make  

# ppm-dpd_patch: 
# 	patch ${WRKDIR}/ppm-dpd/src/dpd_module_random.f90 patches/dpd-a
# 	patch ${WRKDIR}/ppm-dpd/src/dpd_cumDist.f90 patches/dpd-b
