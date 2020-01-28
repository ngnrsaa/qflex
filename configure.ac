#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([qFlex], [0.1], [salvatore.mandra@nasa.gov])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_PREFIX_DEFAULT([$HOME/local])

AC_LANG(C++)
AC_PROG_CXX

# Default CXXFLAGS.
CXXFLAGS="-O3  -std=c++17  -march=native -Wall -Wpedantic -Wno-unused-command-line-argument "$CXXFLAGS

# Default LIBS
LIBS="-lpthread "$LIBS

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/src/Makefile
                 tests/python/Makefile])

# Skip all checks
AC_ARG_ENABLE([all_checks], 
  AS_HELP_STRING([--disable-all_checks],[Disable all the checks.]), 
  [case "${enableval}" in
      yes) all_checks=true  ;;
       no) all_checks=false ;; 
        *) AC_MSG_ERROR([bad value ${enableval} for --disable-all_checks])
   esac], [all_checks=true]
)

# Skip error on warning
AC_ARG_ENABLE([error_on_warnings], 
  AS_HELP_STRING([--disable-error_on_warnings],[Disable error on warnings.]), 
  [case "${enableval}" in
      yes) error_on_warnings=true  ;;
       no) error_on_warnings=false ;; 
        *) AC_MSG_ERROR([bad value ${enableval} for --disable-error_on_warnings])
   esac], [error_on_warnings=true]
)

# Check for additional Mac OSX dependencies
if test $(uname) == "Darwin"; then
  AC_ARG_ENABLE([brew], 
    AS_HELP_STRING([--disable-brew],[Disable Brew.]), 
    [case "${enableval}" in
        yes) brew=true  ;;
         no) brew=false ;; 
          *) AC_MSG_ERROR([bad value ${enableval} for --disable-brew])
     esac], [brew=true]
  )
else
  brew=false
fi

# Skip python checks
AC_ARG_ENABLE([python_checks], 
  AS_HELP_STRING([--disable-python_checks],[Disable Python checks.]), 
  [case "${enableval}" in
      yes) python_checks=true  ;;
       no) python_checks=false ;; 
        *) AC_MSG_ERROR([bad value ${enableval} for --disable-python_checks])
   esac], [python_checks=true]
)

# Add extra param for cirq
AC_ARG_ENABLE([python_tests], 
  AS_HELP_STRING([--disable-python_tests],[Disable Python tests.]), 
  [case "${enableval}" in
      yes) python_tests=true  ;;
       no) python_tests=false ;; 
        *) AC_MSG_ERROR([bad value ${enableval} for --disable-python_tests])
   esac], [python_tests=true]
)

# Add extra param for pybind11
AC_ARG_ENABLE([pybind11], 
  AS_HELP_STRING([--disable-pybind11],[Disable pybind11 extension of qFlex.]), 
  [case "${enableval}" in
      yes) pybind11=true  ;;
       no) pybind11=false ;; 
        *) AC_MSG_ERROR([bad value ${enableval} for --disable-python_tests])
   esac], [pybind11=true]
)

# Add extra param for OpenMP
AC_ARG_ENABLE([openmp], 
  AS_HELP_STRING([--enable-openmp],[(Experimental) Enable OpenMP.]), 
  [case "${enableval}" in
      yes) openmp=true  ;;
       no) openmp=false ;; 
        *) AC_MSG_ERROR([bad value ${enableval} for --enable-openmp])
   esac], [openmp=false]
)

# Check for BREW.
if test ${all_checks} != "false" && ${brew} != "false"; then
  AC_CHECK_PROG(BREW_CHECK,brew,yes)
  AS_IF([test x"$BREW_CHECK" != x"yes"], [AC_MSG_ERROR([Please install brew before configuring.])])

  # Check if dependencies are installed
  for pkg in coreutils openblas autoconf automake libomp pybind11; do
    AC_MSG_CHECKING([for ${pkg}])
    if brew ls --versions ${pkg} >/dev/null; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([Please run 'brew install ${pkg}' to install qFlex.])
    fi
  done

  # Update pre-processor flags
  CPPFLAGS=${CPPFLAGS}" -I$(brew --prefix)/include/ -L$(brew --prefix)/lib/"
fi

if test ${all_checks} != "false"; then

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_AWK
AC_PROG_MAKE_SET

# Check for functions.
AC_CHECK_FUNCS([pow], [], [AC_MSG_ERROR([Function 'pow' is missing.])])

# Check for inline.
AC_C_INLINE

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
  AC_CHECK_LIB([openblas], [cblas_cgemm], [
    AC_CHECK_HEADER([cblas.h], [], [
      AC_MSG_ERROR([OpenBLAS installed but cblas.h not found.])
    ])
  ], [
    AC_MSG_ERROR([Please install OpenBLAS before configuring.])
  ])
  
)

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

AS_IF([test ${openmp} == "true"], 
  AC_CHECK_HEADER([omp.h], [], [AC_MSG_ERROR([Please install omp.h before configuring.])])
  AC_MSG_WARN([OpenMP is still experimental with potential memory leak.])
)

# Check for python. If --disable-pybind11, no reasons to activate python
REQ_PYTHON_VERSION="3.5"
if test ${pybind11} == "true" && ${python_checks} == "true"; then
  AM_PATH_PYTHON([${REQ_PYTHON_VERSION}],, 
    pybind11=false
    AC_MSG_WARN([Please install Python >= ${REQ_PYTHON_VERSION} to enable Python porting of qFlex.]))
fi

python_check_module() {

  # Get module
  module=${1}

  # Get version
  version=${2}

  # Get disable
  disable=${3}

  AC_MSG_CHECKING([for ${module}])
  if ${PYTHON} -c "
import ${module}

try:
  from packaging import version
  from sys import exit
except:
  exit(0)

assert(version.parse(${module}.__version__) >= version.parse('${version}'))
" 2>/dev/null; then
    AC_MSG_RESULT([yes])
    echo "true"
  else
    AC_MSG_RESULT([no])
    AC_MSG_WARN([Please install ${module} >= '${version}' to enable Python porting of qFlex.])
    echo "false"
  fi

}

# Check for Python version
CIRQ_REQ_VERSION="0.6"
CIRQ_ERR_FLAG=pybind11

PYBIND11_REQ_VERSION="2.2.4"
PYBIND11_ERR_FLAG=pybind11

PYTEST_REQ_VERSION="4.5"
PYTEST_ERR_FLAG=python_tests

# Check python modules. If --disable-pybind11, no reasons to activate python
if test ${pybind11} == "true" && ${python_checks} == "true"; then

  if test $(python_check_module cirq     ${CIRQ_REQ_VERSION}     ${CIRQ_ERR_FLAG})    == "false"; then
    pybind11="false"
  fi

  if test $(python_check_module pybind11 ${PYBIND11_REQ_VERSION} ${PYBIND11_ERR_FLAG}) == "false"; then
    pybind11="false"
  fi

  if test $(python_check_module pytest   ${PYTEST_REQ_VERSION}   ${PYTEST_ERR_FLAG})   == "false"; then
    pybind11="false"
  fi

fi

# If --disable-pybind11, switch off everything
AS_IF([test ${pybind11} == "false"],
  python_tests=false
) 

fi

AC_SUBST(openmp, ${openmp})
AC_SUBST(pybind11, ${pybind11})
AC_SUBST(python_tests, ${python_tests})
AC_SUBST(error_on_warnings, ${error_on_warnings})

AC_OUTPUT
