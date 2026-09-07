# ###############################################################################
# i-jedi - Interface to JEDI
# ###############################################################################

cmake_minimum_required(VERSION 3.23 FATAL_ERROR)

project(ijedi VERSION 0.0.1 LANGUAGES C CXX Fortran)

option(OPENMP "Build i-jedi with OpenMP support" ON)

find_package(ecbuild 3.5 REQUIRED)
include(ecbuild_system NO_POLICY_SCOPE)
ecbuild_declare_project()

include(GNUInstallDirs)

# ###############################################################################
# Dependencies
# ###############################################################################
find_package(NetCDF REQUIRED COMPONENTS C CXX Fortran)

# find_package( gsl-lite REQUIRED HINTS $ENV{gsl_list_DIR})
find_package(eckit 1.24.4 REQUIRED)
find_package(fckit 0.11.0 REQUIRED)
find_package(atlas 0.36.0 REQUIRED)


#find_package( vader  1.7.0    REQUIRED)
find_package( oops 1.10.0     REQUIRED)
find_package( mist            REQUIRED) # no release number
find_package( saber  1.10.0   REQUIRED)
find_package( ioda   2.9.0    REQUIRED)
find_package( ufo    1.10.0   REQUIRED)
if(OPENMP)
  find_package(OpenMP)
endif()

# ###############################################################################
# Model Specific Dependencies
# ###############################################################################

# FMS for FV3 and MOM6
find_package(FMS 2023.04 COMPONENTS R4 R8)
add_library(fms ALIAS FMS::fms_r8)

# ###############################################################################
# Sources
# ###############################################################################
add_subdirectory(src)

# ###############################################################################
# CTest Directory
# ###############################################################################
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# Build Doxygen documentation
if(ENABLE_DOCS)
  add_subdirectory(docs)
endif()

if(ECBUILD_INSTALL_FORTRAN_MODULES)
  install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/${CMAKE_CFG_INTDIR} DESTINATION ${INSTALL_INCLUDE_DIR})
endif()

# ###############################################################################
# Finalise configuration
# ###############################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project(NAME ${PROJECT_NAME})

# print the summary of the configuration
ecbuild_print_summary()
