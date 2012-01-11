#! /bin/bash
autoheader configure.ac
aclocal -I m4
autoconf 
automake --foreign --add-missing --copy 
