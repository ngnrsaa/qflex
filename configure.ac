#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([qFlex], [1.0], [salvatore.mandra@nasa.gov])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_PREFIX_DEFAULT([$HOME/local])

# Default CXXFLAGS.
: ${CXXFLAGS="-O3  -std=c++17  -march=native"}

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.
AC_CHECK_LIB([gsl], [cblas_cgemm])
AC_CHECK_LIB([gslcblas], [cblas_cgemm])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/Makefile])

AC_OUTPUT
