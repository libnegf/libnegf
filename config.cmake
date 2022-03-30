#######################################################################################
# libNEGF General configuration file
#######################################################################################

# NOTE:
# Set environmnent variables FC and CC so cmake detects the correct compiler 
# module load should do this if properly set

# MPI option
option(WITH_MPI "Whether MPI-parallelised library should be built" FALSE)

# OPEN MP option
option(WITH_OMP "Whether OMP should be used" TRUE)

# Note, this option does not work yet, so leave it to false
option(WITH_INELASTIC "Whether to build with inelastic scattering" FALSE)

# shared library .so build 
option(BUILD_SHARED_LIBS "Whether the library should be shared" FALSE)

# Whether include modules should be added to installation
option(INSTALL_INCLUDE_FILES "Whether module files and headers should be installed" TRUE)

# Whether the tests should be compiled
option(BUILD_TESTING "Whether the tests should be built" TRUE)

# Whether mpifx should be downloaded 
option(FORCE_MPIFX_DOWNLOAD
  "Force mpifx download from repository (do not search for installed package) " FALSE)

# include directory on installation folder
set(INSTALL_INCLUDE_DIR "libnegf" CACHE PATH
  "Installation directory for C and C++ header files (within standard include folder)")

# fortran modules directory on installation folder
set(INSTALL_MOD_DIR "${INSTALL_INCLUDE_DIR}/modfiles" CACHE PATH
  "Installation directory for Fortran module files (within standard include folder)")

#######################################################################################
# Lapack library, manual setting example
#######################################################################################
#set(MKL_LIBDIR "/usr/pack/intel_mkl-11.1-ma/mkl/lib/intel64" CACHE STRING
#	"Directory where to look for mkl libraries"	)

#set(BLAS_LIBRARY_DIR "${MKL_LIBDIR};/usr/lib/x86_64-linux-gnu"
#	CACHE STRING "Paths to scan when looking for the LAPACK libraries")

#set(BLAS_LIBRARY "mkl_gnu_thread;mkl_core;-lgomp" CACHE STRING "BLAS libraries")

#set(LAPACK_LIBRARY_DIR "${MKL_LIBDIR}"
#	CACHE STRING "Paths to scan when looking for the LAPACK libraries")

#set(LAPACK_LIBRARY "mkl_intel_lp64" CACHE STRING "LAPACK libraries")



