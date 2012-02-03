#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.62)
AC_INIT(metis, 4.1.1)
AC_CONFIG_SRCDIR([src/balance.c])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE(metis, 4.1.1)

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([malloc.h stdlib.h string.h strings.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
