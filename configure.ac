#
AC_INIT([OpenEMTP], [2.0], [],[oemtp],[http://usp.br])

AC_LANG(C++)
#libboost1.46
#libboost-thread1.46
#PKG_CHECK_MODULES([DEPS], [boost >= 1.38 libboost-thread >= 1.38 libc >= 6.0 libncurses >= 5.7 libcomedi >= 0.8])

PKG_CHECK_MODULES([DEPS], [gsl])
#AC_SUBST(DEPS_CFLAGS)
#AC_SUBST(DEPS_LIBS)
CPPFLAGS="$CPPFLAGS $DEPS_CFLAGS" 
LIBS="$DEPS_LIBS $LIBS -lgsl -lgslcblas -lm -lgfortran" 

AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([1.10 -Wall no-define])

AC_CONFIG_HEADERS([config.h])
AC_PROG_CXX
AC_PROG_F77(gfortran)

FFLAGS="-g -O0"
CXXFLAGS="-O0 -g3"
#FFLAGS="-O3"
#CXXFLAGS="-O3"
#LT_INIT([disable-static])
LT_INIT

AC_SUBST([OEMTP_SO_VERSION], [1:2:0])
AC_SUBST([OEMTP_API_VERSION], [2.0])

AC_ARG_WITH([openrelay],
[AS_HELP_STRING([--with-openrelay],[OpenRelay Support])],
[],[with_openrelay=no])
AS_IF(	[test "x$with_openrelay" != xno],
	[
		#Set up OpenRelay
		CPPFLAGS="-DOPENRELAY `pkg-config openrelay-1.0 --cflags --libs` -lboost -lpthread  -lboost_date_time-mt -lboost_thread-mt $CPPFLAGS"
	]
	)
	
AC_ARG_WITH([debug],
[AS_HELP_STRING([--with-debug],[Debug])],
[],[with_debug=no])
AS_IF(	[test "x$with_debug" != xno],
	[
		#Set up OpenRelay
		CPPFLAGS="-DDEBUG $CPPFLAGS"
	]
	)

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi
AM_CONDITIONAL([HAVE_DOXYGEN], 
[test -n "$DOXYGEN"])AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([doc/Doxyfile])])
#doc/Makefile
#oemt-${OEMTP_API_VERSION}.pc:oemtp.pc.in]
#		
AC_CONFIG_FILES([Makefile
		src/Makefile
		doc/Makefile
		oemtp-${OEMTP_API_VERSION}.pc:oemtp.pc.in
		])
AC_OUTPUT
