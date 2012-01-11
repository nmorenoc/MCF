# AM_PATH_FFTW([ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]])
# ----------------------------------------------------------
# set up for FFTW
#
AC_DEFUN([AM_PATH_FFTW],
[
AC_PROVIDE([AM_PATH_FFTW])

AC_ARG_WITH(fftw,
 AC_HELP_STRING([--with-fftw=PFX], [Prefix where FFTW has been installed] ),
 [
    test "$withval" = no && AC_MSG_ERROR([fftw is a required package])
    test "$withval" = yes || fftw_prefix="$withval" 
    with_fftw=yes ],
 [ with_fftw=yes ] )

if test $with_fftw = yes ; then
#user override
AS_IF([test "x$FFTW_LIBS" != x && test "x$FFTW_CXXFLAGS" != x ],
[
  have_fftw=yes
],
[
saved_LIBS="$LIBS"
saved_CXXFLAGS="$CXXFLAGS"
FFTW_LIBS=""
FFTW_CXXFLAGS=""

if test x$fftw_prefix != x; then
	# very likely the majority of cases, we will have been configured with:
	# --with-fftw=/some/thing
	#

 # should be ac_FFTW_CXXFLAGS="-I$FFTW_prefix/include"
 #  
 ac_FFTW_CXXFLAGS="-I$fftw_prefix/include"
 #
 # Similarly for fftw, the uninstalled library position is simply in
 # $fftw_prefix, but the installed is in the standard prefixed subdirectory.
 #
 # SGI compiler CC (CXX=CC) needs -lm to link maths library, but 
 # GCC c++ does not.
 #
 ac_FFTW_LDOPTS="-L$fftw_prefix/lib"
else
 # the compiler looks in the "standard" places for FFTW.  In real life,
 # it would be quite possible that FFTW would not be installed in
 # /usr/include, /usr/lib etc. so the defaults will not usually find
 # the right dependencies.
 ac_FFTW_CXXFLAGS=""
 ac_FFTW_LDOPTS=""
fi #dnl test fftw_prefix

fftwname="fftw3"
rfftwname="fft3f"

AC_MSG_CHECKING([for fftw_print_max_memory_usage in $fftwname])

	LIBS="$ac_FFTW_LDOPTS $saved_LIBS -l$rfftwname -l$fftwname"
	CXXFLAGS="$ac_FFTW_CXXFLAGS $saved_CXXFLAGS"
	#
	# AC_TRY_LINK uses the c compiler (set by AC_LANG), so we will
	# temporarily reassign $CC to the c++ compiler.
 	#
	AC_LANG_PUSH(C++)
	AC_TRY_LINK([#include <$fftwname.h>] ,[  fftw_print_max_memory_usage();  ], have_fftw=yes, have_fftw=no)
	if test x$have_fftw=xyes; then
	   AC_TRY_LINK(
[#include <$fftwname.h>] ,[  
   fftw_real *fftwp = 0;
   float *fftp = 0;
   fftp = fftwp;
          ], 
           have_fftw=yes, have_fftw=no)
	fi
	AC_MSG_RESULT($have_fftw)

if test $have_fftw = no; then

  fftwname="sfftw"
  rfftwname="srfftw"
  AC_MSG_CHECKING([for fftw_print_max_memory_usage in $fftwname])
  LIBS="$ac_FFTW_LDOPTS $saved_LIBS -l$rfftwname -l$fftwname"
  CXXFLAGS="$ac_FFTW_CXXFLAGS $saved_CXXFLAGS"
  AC_TRY_LINK([#include <$fftwname.h>] ,[  fftw_print_max_memory_usage();  ], have_fftw=yes, have_fftw=no)
  if test x$have_fftw=xyes; then
    AC_TRY_LINK(
      [#include <$fftwname.h>] ,[
   fftw_real *fftwp = 0;
   float *fftp = 0;
   fftp = fftwp;
      ],
      have_fftw=yes, have_fftw=no)
    fi
  AC_MSG_RESULT($have_fftw)
fi

 AC_LANG_POP(C++)
 LIBS="$saved_LIBS"
 CXXFLAGS="$saved_CXXFLAGS"

]) #dnl user override

AS_IF([test x$have_fftw = xyes],
  [
    test "x$FFTW_CXXFLAGS" = x && FFTW_CXXFLAGS="$ac_FFTW_CXXFLAGS"
    test "x$FFTW_LIBS" = x && FFTW_LIBS="$ac_FFTW_LDOPTS -l$rfftwname -l$fftwname"
    ifelse([$1], , :, [$1])
  ],
  [
    AC_MSG_ERROR([If fftw exist on you system, are you sure you are using the 
    fftw libraries that was configured with --enable-float?])
    ifelse([$2], , :, [$2])
  ])

fi # --with-fftw

AC_SUBST(FFTW_CXXFLAGS)
AC_SUBST(FFTW_LIBS)

])
