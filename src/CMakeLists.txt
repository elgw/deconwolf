# This file is to create a library out of some of the functions
# in dw.
# $ mkdir build
# $ cd build
# $ cmake ..
# or, to change the root of the installation, use something like:
# $ cmake .. -DCMAKE_INSTALL_PREFIX=/usr/
# $ sudo make install
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

option (ENABLE_GPU "Enable GPU/OpenCL" ON)
option (PYDECONWOLF "Build Python module" ON)

add_library( deconwolf SHARED
  dw_align_dots.c
  dw_background.c
  dw_align_dots.c
  dw_dots.c
  dw_imshift.c
  dw_maxproj.c
  dw_png.c
  dw_psf.c
  dw_psf_sted.c
  dw_tiff_merge.c
  dw_util.c
  dw_bwpsf.c
  bw_gsl.c
  li.c
  fft.c
  fim.c
  fim_tiff.c
  ftab.c
  gmlfit.c
  method_identity.c
  method_rl.c
  method_shb.c
  npio.c
  qalign.c
  quickselect.c
  sparse_preprocess.c
  sparse_preprocess_cli.c
  tiling.c
  dw.c
  lanczos.c
)

set(DECONWOLF_PUBLIC_HEADERS
  "dw_util.h"
  "fim.h"
  "fim_tiff.h"
  "fft.h"
  "ftab.h"
  "dw_version.h"
  "npio.h"
  "dw_bwpsf.h"
  "dw.h"
  "kdtree/include/kdtree.h")

set(CMAKE_BUILD_WITH_INSTALL_RPATH 1)

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
if(NOT WIN32)
  target_link_libraries(deconwolf FFTW::FloatOpenMP)
endif()

#
# LibPNG
#
find_package(PNG)
target_link_libraries(deconwolf PNG::PNG)


#
# kdtree
#
add_subdirectory("kdtree/")
target_include_directories(deconwolf PUBLIC "kdtree/include/")
target_link_directories(deconwolf PUBLIC "kdtree/")
target_link_libraries(deconwolf kdtree)

#
# trafo
#
add_subdirectory("trafo")
target_include_directories(deconwolf PUBLIC "trafo/src/")
target_link_directories(deconwolf PUBLIC "trafo/")
target_link_libraries(deconwolf trafo)

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


#
# OpenCL
#
# https://cmake.org/cmake/help/latest/module/FindOpenCL.html

if (ENABLE_GPU)
  find_package(OpenCL)
  if(OpenCL_FOUND)
    message("Will use OpenCL")
    target_link_libraries(deconwolf OpenCL::OpenCL)
    add_definitions(-DOPENCL)
    add_definitions(-DVKFFT_BACKEND=3)
    target_include_directories(deconwolf PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/VkFFT/vkFFT/)
    target_sources(deconwolf PRIVATE cl_util.c)
    target_sources(deconwolf PRIVATE method_shb_cl.c)
    target_sources(deconwolf PRIVATE method_shb_cl2.c)
  else()
    message("No OpenCL :(")
  endif()
endif()

if (PYDECONWOLF)
  add_subdirectory(pywrapper)
endif()



include("GNUInstallDirs")

set_target_properties( deconwolf PROPERTIES
  PUBLIC_HEADER "${DECONWOLF_PUBLIC_HEADERS}")

INSTALL(TARGETS deconwolf
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_FULL_INCLUDEDIR}/deconwolf"
  LIBRARY DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
)
