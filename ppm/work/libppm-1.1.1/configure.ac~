#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(ppm, 1.1.1)
AC_CONFIG_SRCDIR([src/ppm_define.h])
AC_CONFIG_HEADER([config.h])
AM_INIT_AUTOMAKE(libppm, 1.1.1)

# Checks for programs.
AC_PROG_FC
AC_LANG(Fortran)

# use mpi in any case (probably fakempi.h)
AC_DEFINE([USE_MPI], [1], [use mpi])
AC_ARG_ENABLE(mpi, [AC_HELP_STRING([--enable-mpi],[compile  with mpi])], have_mpi=$enableval, have_mpi=no)
if test "$have_mpi" = "yes"; then
   ACX_MPI(,
	   [AC_MSG_ERROR([cannot find working mpif90, please set MPIFC variable])])
   FC="$MPIFC"
   LIBS="$MPILIBS $LIBS"
fi

AM_CONDITIONAL([HAVE_MPI], [test x$have_mpi = xyes])


AC_PROG_RANLIB
AC_CHECK_PROG(MAKEDEPF90, [makedepf90], [makedepf90], [not found], [$PATH])
AC_ARG_VAR(MAKEDEPF90, [F90 dependency tracking])

# HACK: should be changed to more aproprite linking of c++ and fortran code
AC_PROG_CXX
AC_LANG_PUSH(C++)
AC_CHECK_LIB(stdc++,main,,AC_MSG_ERROR(Cannot find lib stdc++))
AC_CHECK_LIB(gcc_s,main,,AC_MSG_WARN(Cannot find lib gcc_s))
AC_LANG_POP

# Checks for libraries.
LIBS="$LIBS $MPILIBS"
AC_CHECK_LIB(vizing, 
	     vizing_coloring,
	     ,
	     [AC_MSG_ERROR([Cannot find lib vizing, please set LDFLAGS])])

AC_CHECK_LIB(metis, 
	     METIS_WPartGraphRecursive,
	     [],
	     [AC_MSG_ERROR([Cannot find lib metis, please set LDFLAGS])])

#AC_LANG_PUSH(C)
#CPPFLAGS=$FCFLAGS
#AC_CHECK_HEADER(fftw3.h, [], [AC_MSG_ERROR([cannot find fftw3.h, please set FCFLAGS])])
#AC_LANG_POP

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

AX_WITH_PROG([FUNIT],[funit],[not])
AM_CONDITIONAL([HAVE_FUNIT], [! test x$FUNIT = xnot])

# Checks for library functions.
AC_CONFIG_FILES([Makefile src/Makefile test/Makefile])
AC_OUTPUT


AC_MSG_NOTICE([================================================])
AC_MSG_NOTICE([FC=$FC])
AC_MSG_NOTICE([FCFLAGS=$FCFLAGS])
AC_MSG_NOTICE([LIBS=$LIBS])
AC_MSG_NOTICE([MAKEDEPF90=$MAKEDEPF90])
if test "$have_mpi" = "yes"; then
   AC_MSG_NOTICE([using mpi])
else 
   AC_MSG_NOTICE([not using mpi])
fi
if test x$FUNIT = xnot; then
   AC_MSG_NOTICE([not using funit])
else 
   AC_MSG_NOTICE([using funit])
   AC_MSG_NOTICE([FUNIT=$FUNIT])
fi
