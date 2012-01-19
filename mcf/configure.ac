#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


#Initialize Autoconf.
#Specify package's name, version number and bug-report email address.
AC_INIT([mcf], [1.0], [xin.bian@aer.mw.tum.de])

#Require a minimum Autoconf version.
AC_PREREQ(2.60)

#A safety check. 
#src/mcf.F90 should be a distributed source file.
AC_CONFIG_SRCDIR([src/mcf.F90])

#Auxiliary scripts such as install-sh and depcomp 
#should be in auxilary folder.
AC_CONFIG_AUX_DIR([build-aux])


#Declear config.h as output header.
AC_CONFIG_HEADERS([config.h])

#Initialize Automake.
#Turn on all Automake warnings and report them as errors.
#This is a foreign package.
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

#Compiler checks.
AC_PROG_FC
AC_LANG(Fortran)
ACX_MPI(,
	[AC_MSG_ERROR([cannot find working mpif90, please set MPIFC variable])])

# compile and link with mpif90
CC="$MPIFC"
FC="$MPIFC"
LIBS="$MPILIBS $LIBS"

# HACK: should be changed to more aproprite linking of c++ and fortran code
AC_PROG_CXX
AC_LANG_PUSH(C++)
AC_CHECK_LIB(stdc++,main,,AC_MSG_ERROR(Cannot find lib stdc++))
AC_CHECK_LIB(gcc_s,main,,AC_MSG_WARN(Cannot find lib gcc_s))
AC_LANG_POP

# for dependencies tracking
AC_CHECK_PROG(MAKEDEPF90, [makedepf90], [makedepf90], [not found], [$PATH])
AC_ARG_VAR(MAKEDEPF90, [F90 dependency tracking])


# ppm library
#ac_ext=f90 

AC_FC_SRCEXT(F90)


# Checks for libraries.

# double precision fftw library
AC_CHECK_LIB(fftw3, 
	     dfftw_execute_dft_c2r,
	     [],
	     [AC_MSG_ERROR([cannot link against fftw library, please set LDFLAGS])])

# single precision fftw library
AC_CHECK_LIB(fftw3f, 
	     sfftw_execute_dft_c2r,
	     [],
	     [AC_MSG_ERROR([cannot link against fftw library, please set LDFLAGS])])

# vizing library
AC_CHECK_LIB(vizing, 
	     vizing_coloring,
	     ,
	     [AC_MSG_ERROR([Cannot find lib vizing, please set LDFLAGS])])
# metis library
AC_CHECK_LIB(metis, 
	     METIS_WPartGraphRecursive,
	     [],
	     [AC_MSG_ERROR([Cannot find lib metis, please set LDFLAGS])])

LIBS="-lppm $LIBS"

AC_LINK_IFELSE([
PROGRAM main
	USE ppm_module_write
	INTEGER :: i
	CALL ppm_write(1, 'test', 'test', i)
END PROGRAM main
],
	[AC_MSG_RESULT([can compile with ppm])],
	[AC_MSG_ERROR([cannot link against ppm library, please set LDFLAGS and FCFLAGS])])


#Declare Makefile and src/Makefile as output files
AC_CONFIG_FILES([Makefile src/Makefile])

#Actually output all declared files.
AC_OUTPUT
