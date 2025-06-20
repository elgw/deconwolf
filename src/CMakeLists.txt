# This file is to create a library out of some of the functions
# in dw.
# mkdir build
# cd build
# cmake .. -DCMAKE_INSTALL_PREFIX="/"
# sudo make install
#
# Using some other install prefix,
# you might want to add the path to the file
# `/etc/ld.so.conf`
# and then run `ldconfig`

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(libdw LANGUAGES C)

set (CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

option (ENABLE_GPU "Enable GPU/OpenCL" OFF)

add_library( deconwolf SHARED
  fim.c
  fim_tiff.c
  ftab.c
  fft.c
  dw_util.c
  gmlfit.c
  lanczos.c
  tiling.c
  quickselect.c
  npio.c
)

set(DECONWOLF_PUBLIC_HEADERS
  "dw_util.h"
  "fim.h"
  "fim_tiff.h"
  "fft.h"
  "ftab.h"
  "dw_version.h"
  "npio.h"
"kdtree/include/kdtree.h")

#
# Gnu Scientific Library
#
find_package(GSL)
target_link_libraries(deconwolf ${GSL_LIBRARIES})

# FFTW3
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../util/findFFTW)
if(UNIX)
  find_package(FFTW REQUIRED FLOAT_LIB FLOAT_OPENMP_LIB)
else()
  find_package(FFTW REQUIRED FLOAT_LIB)
endif()
target_link_libraries(deconwolf FFTW::Double)
target_link_libraries(deconwolf FFTW::Float)


#
# kdtree
#
add_subdirectory("kdtree/")
target_include_directories(deconwolf PUBLIC "kdtree/include/")
target_link_directories(deconwolf PUBLIC "kdtree/")
target_link_libraries(deconwolf kdtree)


#
# OpenMP
#
#  https://cmake.org/cmake/help/v3.9/module/FindOpenMP.html
find_package(OpenMP)
if(OpenMP_C_FOUND)
  message("Found OpenMP")
  target_link_libraries(deconwolf OpenMP::OpenMP_C)
else()
  message("OpenMP not found automatically, trying manual")

  set(
    CMAKE_C_FLAGS
    "${CMAKE_C_FLAGS} -fopenmp=omp /openmp"
  )
  SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fopenmp=omp -lomp")

endif()

#
# TIFF
#
find_package(TIFF)
target_link_libraries(deconwolf ${TIFF_LIBRARIES})

include("GNUInstallDirs")

set_target_properties( deconwolf PROPERTIES
PUBLIC_HEADER "${DECONWOLF_PUBLIC_HEADERS}")

INSTALL(TARGETS deconwolf
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_FULL_INCLUDEDIR}/deconwolf"
  LIBRARY DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
  )
