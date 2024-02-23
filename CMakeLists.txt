# This is a file to be interpreted by cmake (https://cmake.org/)
#
# Typically it is run to generate a makefile (for make) like this:
#
# mkdir build
# cd build
# cmake ..
# make
#
# Variations:
# - To Specify the compiler:
#   CC=clang cmake ..
# - To see command lines:
#   make VERBOSE=1
#
# Useful references:
# Good tutorial:
# https://cliutils.gitlab.io/modern-cmake/
# Complex example:
# https://github.com/fluent/fluent-bit/blob/master/CMakeLists.txt


cmake_minimum_required(VERSION 3.9)


if(WIN32)
    # Set MSVC runtime to MultiThreaded (/MT)
    cmake_policy(SET CMP0091 NEW)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

project(deconwolf VERSION 0.3.6
  DESCRIPTION "Deconvolution of wide field images"
  LANGUAGES C)
set (CMAKE_C_STANDARD 99)

include(CheckIPOSupported)

if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
  message("This is a ${CMAKE_SYSTEM_NAME} system")
elseif(${CMAKE_SYSTEM_NAME} MATCHES Windows)
  message("This is a Windows System")
endif()

# UNIX, WIN32, WINRT, CYGWIN, APPLE are environment variables as flags set by default system
if(UNIX)
  message("This is a ${CMAKE_SYSTEM_NAME} system")
  add_definitions(-DUNIX)
elseif(WIN32)
  message("This is a Windows System")
  add_definitions(-DWINDOWS)
endif()


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# here is where we set the language(deconwolf C)
add_executable(dw
  src/deconwolf.c
  src/dw.c
  src/dw_dots.c
  src/dw_imshift.c
  src/dw_maxproj.c
  src/dw_psf.c
  src/dw_psf_sted.c
  src/dw_tiff_merge.c
  src/dw_util.c
  src/fim.c
  src/fim_tiff.c
  src/ftab.c
  src/fwhm.c
  src/method_identity.c
  src/method_rl.c
  src/method_shb.c
  src/tiling.c
  src/fft.c)

add_executable(dw_bw
  src/bw_gsl.c
  src/dw_bwpsf.c
  src/dw_util.c
  src/fft.c
  src/fim.c
  src/fim_tiff.c
  src/li.c
  src/ftab.c
  src/lanczos.c)

#
# FFTW3
#

target_link_libraries(dw fftw3 fftw3f fftw3f_omp)
target_link_libraries(dw_bw fftw3 fftw3f fftw3f_omp)

#
# Math library, if needed
#

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
  target_link_libraries(dw ${MATH_LIBRARY})
  target_link_libraries(dw_bw ${MATH_LIBRARY})
endif()

#
# TIFF
#

target_link_libraries(dw tiff)
target_link_libraries(dw_bw tiff)

#
# Gnu Scientific Library
#

target_link_libraries(dw gsl)
target_link_libraries(dw_bw gsl)

#
# OpenCL
#
# https://cmake.org/cmake/help/latest/module/FindOpenCL.html

find_package(OpenCL)
if(OpenCL_FOUND)
  message("Will use OpenCL")
  target_link_libraries(dw OpenCL)
  add_definitions(-DOPENCL)
  add_definitions(-DVKFFT_BACKEND=3)
  target_include_directories(dw PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/VkFFT/vkFFT/)
  target_sources(dw PRIVATE src/cl_util.c)
  target_sources(dw PRIVATE src/method_shb_cl.c)
  target_sources(dw PRIVATE src/method_shb_cl2.c)
else()
  message("No OpenCL :(")
endif()



#
# OpenMP
#
#  https://cmake.org/cmake/help/v3.9/module/FindOpenMP.html
if(WIN32)
  message("Enabling OpenMP since it should be available")
  add_definitions(-D_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp=libomp")
else()
find_package(OpenMP)
if(OpenMP_C_FOUND)
  message("Found OpenMP")
  target_link_libraries(dw OpenMP::OpenMP_C)
  target_link_libraries(dw_bw OpenMP::OpenMP_C)
endif()
endif()


#include(CheckCCompilerFlag)
#check_c_compiler_flag(-flto=auto CP_ACCEPT_FLTOAUTO)

# Enable -flto (GCC) and similar if available, see
# https://cmake.org/cmake/help/latest/policy/CMP0069.html
check_ipo_supported(RESULT result)
if(result)
  set_property(DIRECTORY PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
