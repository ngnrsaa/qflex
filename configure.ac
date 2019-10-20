#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([qFlex], [1.0], [salvatore.mandra@nasa.gov])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_PREFIX_DEFAULT([$HOME/local])

# Default CXXFLAGS.
: ${CXXFLAGS="-O3  -std=c++17  -march=native"}

# Checks for programs.
AC_LANG(C++)
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Check for GIT.
AC_CHECK_PROG(GIT_CHECK,git,yes)
AS_IF([test x"$GIT_CHECK" != x"yes"], [AC_MSG_ERROR([Please install git before configuring.])])

# Checks for headers.
AC_CHECK_HEADER([algorithm], [], [AC_MSG_ERROR([Please install algorithm before configuring.])])
AC_CHECK_HEADER([cassert], [], [AC_MSG_ERROR([Please install cassert before configuring.])])
AC_CHECK_HEADER([chrono], [], [AC_MSG_ERROR([Please install chrono before configuring.])])
AC_CHECK_HEADER([cmath], [], [AC_MSG_ERROR([Please install cmath before configuring.])])
AC_CHECK_HEADER([complex], [], [AC_MSG_ERROR([Please install complex before configuring.])])
AC_CHECK_HEADER([ctime], [], [AC_MSG_ERROR([Please install ctime before configuring.])])
AC_CHECK_HEADER([fstream], [], [AC_MSG_ERROR([Please install fstream before configuring.])])
AC_CHECK_HEADER([iostream], [], [AC_MSG_ERROR([Please install iostream before configuring.])])
AC_CHECK_HEADER([iterator], [], [AC_MSG_ERROR([Please install iterator before configuring.])])
AC_CHECK_HEADER([list], [], [AC_MSG_ERROR([Please install list before configuring.])])
AC_CHECK_HEADER([memory], [], [AC_MSG_ERROR([Please install memory before configuring.])])
AC_CHECK_HEADER([omp.h], [], [AC_MSG_ERROR([Please install omp.h before configuring.])])
AC_CHECK_HEADER([regex], [], [AC_MSG_ERROR([Please install regex before configuring.])])
AC_CHECK_HEADER([sstream], [], [AC_MSG_ERROR([Please install sstream before configuring.])])
AC_CHECK_HEADER([string], [], [AC_MSG_ERROR([Please install string before configuring.])])
AC_CHECK_HEADER([unordered_map], [], [AC_MSG_ERROR([Please install unordered_map before configuring.])])
AC_CHECK_HEADER([unordered_set], [], [AC_MSG_ERROR([Please install unordered_set before configuring.])])
AC_CHECK_HEADER([vector], [], [AC_MSG_ERROR([Please install vector before configuring.])])

# Checks for libraries.
AS_IF([test x${CXX} == x"icpc"],
  AC_CHECK_HEADER([mkl.h], [], [AC_MSG_ERROR([Please install mkl.h before configuring.])])
  ,
  AC_CHECK_HEADER([gsl/gsl_cblas.h], [], [AC_MSG_ERROR([Please install gsl/gsl_cblas.h before configuring.])])
  AC_CHECK_LIB([gslcblas], [cblas_cgemm], [], [AC_MSG_ERROR([Please install gslcblas before configuring.])])
)

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/src/Makefile
                 tests/python/Makefile])

# Add extra param for cirq
AC_ARG_ENABLE([cirq], 
  AS_HELP_STRING([--disable-cirq],[Disable cirq dependant features.]), 
  [case "${enableval}" in
      yes) cirq=true  ;;
       no) cirq=false ;; 
        *) AC_MSG_ERROR([bad value ${enableval} for --disable-cirq])
   esac], [cirq=true]
)

# Add extra param for pybind11
AC_ARG_ENABLE([pybind11], 
  AS_HELP_STRING([--disable-pybind11],[Disable pybind11 extension of qFlex.]), 
  [case "${enableval}" in
      yes) pybind11=true  ;;
       no) pybind11=false ;; 
        *) AC_MSG_ERROR([bad value ${enableval} for --disable-cirq])
   esac], [pybind11=true]
)

# Check for Python version
REQ_PYTHON_VERSION="3.5"
AM_PATH_PYTHON([${REQ_PYTHON_VERSION}],, 
  AS_IF([test ${cirq} == "true" || ${pybind11} == "true"], AC_MSG_ERROR([Please install Python >= ${REQ_PYTHON_VERSION} before configuring.]))
)

# Check that cirq is installed
REQ_CIRQ_VERSION="0.5"
AS_IF([test ${cirq} == "true"],
  AC_MSG_CHECKING([for cirq])
  AS_IF([${PYTHON} -c "import cirq; import sys; from packaging import version; sys.exit(not version.parse(cirq.__version__) >= version.parse('${REQ_CIRQ_VERSION}'))" 2>/dev/null],
    AC_MSG_RESULT([yes]), 
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([Please install cirq >= '${REQ_CIRQ_VERSION}' before configuring (or --disable-cirq).]))
)

AC_SUBST(cirq, ${cirq})

# Check that pybind11 is installed
REQ_PYBIND11_VERSION="2.2.4"
AS_IF([test ${pybind11} == "true"],
  AC_MSG_CHECKING([for pybind11])
  AS_IF([${PYTHON} -c "import pybind11; import sys; from packaging import version; sys.exit(not version.parse(pybind11.__version__) >= version.parse('${REQ_PYBIND11_VERSION}'))" 2>/dev/null],
    AC_MSG_RESULT([yes]), 
    AC_MSG_RESULT([no])
    AC_MSG_ERROR([Please install pybind11 >= '${REQ_PYBIND11_VERSION}' before configuring (or --disable-pybind11).]))
)

AC_SUBST(pybind11, ${pybind11})

# Check that pytest is installed
REQ_PYTEST_VERSION="4.5"
AS_IF([test ${pybind11} == "true"],
  AS_IF([test ${cirq} == "true"],
    AC_MSG_CHECKING([for pytest])
    AS_IF([${PYTHON} -c "import pytest; import sys; from packaging import version; sys.exit(not version.parse(pytest.__version__) >= version.parse('${REQ_PYTEST_VERSION}'))" 2>/dev/null],
      AC_MSG_RESULT([yes])
      pytest=true, 
      AC_MSG_RESULT([no])
      pytest=false
      AC_MSG_WARN([Please install pytest >= '${REQ_PYTEST_VERSION}' to run python tests.])),
    pytest=false
  ),
  pytest=false
)

AC_SUBST(pytest, ${pytest})

AC_OUTPUT
