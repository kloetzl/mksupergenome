AC_INIT([mksupergenome], [0.1])
AM_INIT_AUTOMAKE([-Wall foreign ])

AC_CONFIG_MACRO_DIR([m4])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_MAKE_SET
AC_PROG_CPP
AC_PROG_RANLIB
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

# Make sure, also the C++ programs are compiled with OpenMP
AC_LANG(C++)
AC_OPENMP

# Execute all tests using C
AC_LANG(C)
AC_OPENMP

AC_CHECK_LIB([m],[cos])


# By default try to build with libdivsufsort.
AC_ARG_WITH([libdivsufsort],
    AS_HELP_STRING([--without-libdivsufsort], [Build without libdivsufsort and use psufsort instead.]))

AS_IF([test "x$with_libdivsufsort" != "xno"],
	[
		# The libdivsufsort header contains some Microsoft extension making
		# compilation fail on certain systems (i.e. OS X). Add the following
		# flag so the build runs smoothly.
		CPPFLAGS="$CPPFLAGS -fms-extensions"
		AC_CHECK_HEADERS([divsufsort.h],[have_libdivsufsort=yes],[have_libdivsufsort=no])
		AC_CHECK_LIB(divsufsort, divsufsort, [], [have_libdivsufsort=no])
	],
	[
		have_libdivsufsort=no
		# psufsort needs C++11
		AX_CXX_COMPILE_STDCXX_11([],[mandatory])
	]
)

AS_IF([test "x$have_libdivsufsort" = "xno"],
	[AS_IF([test "x$with_libdivsufsort" != "xno"],
		[AC_MSG_ERROR([Configuration for libdivsufsort failed. Either install libdivsufsort, or use our replacement, psufsort, instead.
	./configure --without-libdivsufsort
The latter may result in longer runtimes.])
	])
])

AM_CONDITIONAL([BUILD_WITH_LIBDIVSUFSORT],[test "x${with_libdivsufsort}" != "xno"])


# Check for various headers including those used by libdivsufsort.
AC_CHECK_HEADERS([limits.h stdlib.h string.h unistd.h stdint.h inttypes.h err.h errno.h fcntl.h])

AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_INT32_T
AC_TYPE_UINT8_T
AC_HEADER_STDBOOL

# Until someone convinces me otherwise, I will deactivate the macros
# AC_FUNC_MALLOC and AC_FUNC_REALLOC. They only check if `malloc(0)` retuns a
# non-null pointer. This breaks the build on systems using uClibc, including
# my laptop.
# As requesting zero bytes is not useful, and implementation-defined behaviour,
# it should be avoided in the first place. Thus I really don't need these checks.

AC_CHECK_FUNCS([floor pow sqrt strdup strerror])
AC_CHECK_FUNCS([strtoul strtod])
AC_CHECK_FUNCS([reallocarray])

AM_CONDITIONAL([HAVE_REALLOCARRAY], [test "x$ac_cv_func_reallocarray" = xyes])

AC_CONFIG_HEADERS([config.h:config.hin])

AC_CONFIG_FILES([
 Makefile
 libs/Makefile
])
AC_OUTPUT

