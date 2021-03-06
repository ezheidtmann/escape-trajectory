# ------------------------------------------------------------------------
# $Revision: 1.33 $
# $Date: 2006/03/23 22:44:56 $
# ------------------------------------------------------------------------
# Programmer(s): Radu Serban and Aaron Collier @ LLNL
# ------------------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see sundials/shared/LICENSE.
# ------------------------------------------------------------------------
# Process this file with autoconf to produce a configure script.
# ------------------------------------------------------------------------

# Initialize configure
AC_INIT(SUNDIALS, 2.2.0, radu@llnl.gov)
LIB_REVISION="0:0:0"

# Display greeting
SUNDIALS_GREETING

# Miscellaneous SUNDIALS initializations
echo "Initialization"
echo "--------------"
echo ""
SUNDIALS_INITIALIZE

# Test enable/disable features
SUNDIALS_ENABLES
SUNDIALS_MPI_ENABLES

# Set C compilation (Required)
echo ""
echo "C Compiler Settings"
echo "-------------------"
echo ""
AC_PROG_CC(${CC} cc gcc)
SUNDIALS_SET_CC

# Set Fortran compilation (Optional)
if test "X${F77_ENABLED}" = "Xyes" || test "X${FCMIX_ENABLED}" = "Xyes" ; then
  echo ""
  echo "Fortran Compiler Settings"
  echo "-------------------------"
  echo ""
  AC_PROG_F77(${F77} f77 g77)
  SUNDIALS_SET_F77
fi

# Set MPI support (Optional)
if test "X${MPI_ENABLED}" = "Xyes"; then
  echo ""
  echo "MPI-C Settings"
  echo "--------------"
  echo ""
  SUNDIALS_SET_MPICC
  if test "X${F77_ENABLED}" = "Xyes"; then
    echo ""
    echo "MPI-Fortran Settings"
    echo "--------------------"
    echo ""
    SUNDIALS_SET_MPIF77
  fi
fi

# Run libtool checks
echo ""
echo "Libtool Settings"
echo "----------------"
echo ""
AC_DISABLE_SHARED
AC_LIBTOOL_WIN32_DLL
AC_PROG_LIBTOOL

# Set examples modules
echo ""
echo "Examples"
echo "--------"
echo ""
SUNDIALS_ENABLE_EXAMPLES

# Substitute variables
AC_SUBST(COMPILER_PREFIX)
AC_SUBST(LINKER_PREFIX)
AC_SUBST(OBJ_EXT)
AC_SUBST(MPICC)
AC_SUBST(MPIF77)
AC_SUBST(MPI_INC_DIR)
AC_SUBST(MPI_LIB_DIR)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_FLAGS)
AC_SUBST(MODULES)
AC_SUBST(NVEC_MODULES)
AC_SUBST(EXAMPLE_MODULES)
AC_SUBST(FCMIX_ENABLED)
AC_SUBST(FLOAT_TYPE)
AC_SUBST(LIBTOOL_CMD)
AC_SUBST(LIBTOOL_DEPS)
AC_SUBST(F77_MANGLE_MACRO1)
AC_SUBST(F77_MANGLE_MACRO2)
AC_SUBST(F77_CASE)
AC_SUBST(F77_UNDERSCORES)
AC_SUBST(PRECISION_LEVEL)
AC_SUBST(GENERIC_MATH_LIB)
AC_SUBST(F77_MPI_COMM_F2C)
AC_SUBST(F77_LNKR)
AC_SUBST(MPIF77_LNKR)
AC_SUBST(LIB_REVISION)

# Set installation path
echo ""
echo "Installation Directories"
echo "------------------------"
echo ""
SUNDIALS_INSTALL_PATH

# Generate makefiles
echo ""
echo "Create Makefiles"
echo "----------------"
echo ""
SUNDIALS_BUILD_MODULES_LIST
AC_CONFIG_FILES([${SUNDIALS_MAKEFILES}])
AC_OUTPUT

# Display summary
SUNDIALS_REPORT
