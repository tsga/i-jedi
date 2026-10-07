cmake_minimum_required(VERSION 3.23 FATAL_ERROR)

project(ijedi-bundle VERSION 0.0.1 LANGUAGES C CXX Fortran)

find_package(ecbuild 3.5 REQUIRED)
include(ecbuild_bundle)

ecbuild_bundle_initialize()

# ===============================================================================
# user options
# ===============================================================================
option(TEST_IODA "download and ioda-data repo needed for ioda ctests" OFF)
option(TEST_UFO "download and ufo-data repo needed for ufo ctests" OFF)
option(TEST_VADER "download and jedi-model-data repo needed for saber and vader ctests" OFF)

# ECMWF libraries
option("BUNDLE_SKIP_ECKIT" "Don't build eckit" "ON") # Skip eckit build unless user passes -DBUNDLE_SKIP_ECKIT=OFF
option("BUNDLE_SKIP_FCKIT" "Don't build fckit" "ON") # Skip fckit build unless user passes -DBUNDLE_SKIP_FCKIT=OFF
option("BUNDLE_SKIP_ATLAS" "Don't build atlas" "ON") # Skip atlas build unless user passes -DBUNDLE_SKIP_ATLAS=OFF

# turn off optional OOPS toy models
option(ENABLE_LORENZ95_MODEL "Build Lorenz95 toy model" OFF)
option(ENABLE_QG_MODEL "Build QG toy model" OFF)

# Options to build JCSDA repositories
option(BUILD_JEDI "Build JEDI" ON)
option(CLONE_MIST "Clone mist with ecbuild" ON)

ecbuild_bundle(PROJECT eckit GIT "https://github.com/ecmwf/eckit.git" TAG 1.16.0)
ecbuild_bundle(PROJECT fckit GIT "https://github.com/ecmwf/fckit.git" TAG 0.9.2)
ecbuild_bundle(PROJECT atlas GIT "https://github.com/ecmwf/atlas.git" TAG 0.35.0)

# ===============================================================================
# optional test data (needed first since other bundles check for data)
# ===============================================================================
if(TEST_IODA)
  ecbuild_bundle(PROJECT ioda-data GIT "https://github.com/JCSDA-internal/ioda-data.git" BRANCH develop UPDATE)
endif()

if(TEST_UFO)
  ecbuild_bundle(PROJECT ufo-data GIT "https://github.com/JCSDA-internal/ufo-data.git" BRANCH develop UPDATE)
endif()

if(TEST_VADER)
  ecbuild_bundle(PROJECT jedi-model-data GIT "https://github.com/JCSDA-internal/jedi-model-data.git" BRANCH develop UPDATE)
endif()

# ===============================================================================
# required repositories

# ===============================================================================
ecbuild_bundle( PROJECT jedicmake GIT    "https://github.com/jcsda/jedi-cmake.git"      UPDATE BRANCH develop RECURSIVE )
if (BUILD_JEDI)
  ecbuild_bundle( PROJECT oops      GIT    "https://github.com/tsga/oops.git"          UPDATE BRANCH landVector )
  ecbuild_bundle( PROJECT vader     GIT    "https://github.com/JCSDA/vader.git"         UPDATE BRANCH develop)
  if (CLONE_MIST)
    ecbuild_bundle( PROJECT mist      GIT    "https://github.com/JCSDA-internal/mist.git" UPDATE BRANCH develop )
  else()
    ecbuild_bundle( PROJECT mist SOURCE "./mist" )
  endif()
  ecbuild_bundle( PROJECT saber     GIT    "https://github.com/JCSDA/saber.git"         UPDATE BRANCH develop )
  ecbuild_bundle( PROJECT ioda      GIT    "https://github.com/JCSDA/ioda.git"          UPDATE BRANCH develop )
  ecbuild_bundle( PROJECT ufo       GIT    "https://github.com/JCSDA/ufo.git"           UPDATE BRANCH develop )
  option(BUILD_FV3JEDI "Build fv3-jedi" OFF)
  if (BUILD_FV3JEDI)
    ecbuild_bundle( PROJECT fv3-jedi-lm    GIT "https://github.com/JCSDA/fv3-jedi-linearmodel.git"   UPDATE BRANCH develop )
    option(ENABLE_FV3_JEDI_DATA "Obtain fv3-jedi test data from fv3-jedi-data repository (vs tarball)" ON)
    ecbuild_bundle( PROJECT fv3-jedi-data  GIT "https://github.com/JCSDA-internal/fv3-jedi-data.git" UPDATE BRANCH develop )
    ecbuild_bundle( PROJECT fv3-jedi  GIT    "https://github.com/JCSDA/fv3-jedi.git"   UPDATE BRANCH develop )
  endif()
endif()
ecbuild_bundle(PROJECT ijedi SOURCE "../")

# ===============================================================================
# required repositories
# ===============================================================================
include(jedicmake/cmake/Functions/git_functions.cmake)

# ===============================================================================
ecbuild_bundle_finalize()
