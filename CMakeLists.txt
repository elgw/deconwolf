#
#  This is a file to be interpreted by cmake (https://cmake.org/)
#  Please note that all dependencies have to be available first.
#  See also INSTALL.md
#
# Typical Usage
#
#   mkdir build
#   cd build
#   cmake ..
#   cmake --build .
#
#  On Windows, use:
#
#   mkdir build
#   cd build
#   cmake `-T ClangCL -A x64` ..
#   cmake --build .
#
# Variations/Options
# - On Windows, use
# - To Specify the compiler:
#   CC=clang cmake ..
# - To see command lines:
#   cmake -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON ..
# - See what is defined with
#   cmake -L . or cmake --LAH .
#
#
# Documentation and tutorials:
# https://cmake.org/cmake/help/latest/
# https://cliutils.gitlab.io/modern-cmake/
#
#
# Install generators
# Debian
# https://cmake.org/cmake/help/latest/cpack_gen/deb.html
# RPM
# https://cmake.org/cmake/help/latest/cpack_gen/rpm.html
# Windows
# https://cmake.org/cmake/help/latest/cpack_gen/nsis.html
# Macos
# https://cmake.org/cmake/help/latest/cpack_gen/bundle.html
#

# TODO:
# Check https://learn.microsoft.com/en-us/cpp/error-messages/compiler-warnings/compiler-warning-level-3-c4996?view=msvc-170

cmake_minimum_required(VERSION 3.9)

project(deconwolf VERSION 0.3.7
  DESCRIPTION "Deconvolution of wide field images"
  LANGUAGES C)

set (CMAKE_C_STANDARD 11)

option (ENABLE_GPU "Enable GPU/OpenCL" ON)
option (ENABLE_NATIVE_OPTIMIZATION "Enable non-portable optimizations" OFF)

# UNIX, WIN32, WINRT, CYGWIN, APPLE are environment variables
# as flags set by default system
if(UNIX)
  add_definitions(-DUNIX)
elseif(WIN32)
  add_definitions(-DWINDOWS)
endif()

#
# Default build type is RELEASE
#
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

#
# Add source files
#
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
# Set MSVC runtime to MultiThreaded (/MT)
#

if(WIN32)
  cmake_policy(SET CMP0091 NEW)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

#
# Link time optimization
#

# https://cmake.org/cmake/help/latest/module/CheckIPOSupported.html
include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
  message("Enabling IPO")
  set_property(TARGET dw PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET dw_bw PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_JOBS 8)
else()
  message(WARNING "IPO is not supported: ${output}")
endif()

#
# Architecture optimizations
#
if(ENABLE_NATIVE_OPTIMIZATION)
include(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
endif()
CHECK_C_COMPILER_FLAG("-mtune=native" COMPILER_SUPPORTS_MTUNE_NATIVE)
if(COMPILER_SUPPORTS_MTUNE_NATIVE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mtune=native")
endif()
endif()

#
# FFTW3
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/util/findFFTW)
find_package(FFTW REQUIRED FLOAT_LIB FLOAT_OPENMP_LIB)
target_link_libraries(dw FFTW::Double)
target_link_libraries(dw FFTW::Float)
target_link_libraries(dw FFTW::FloatOpenMP)


# There is no _threads version or _omp version
# on windows, multicore support should be baked
# into the basic library.
if(NOT WIN32)
  target_link_libraries(dw FFTW::FloatOpenMP)
endif()

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
find_package(TIFF)
target_link_libraries(dw ${TIFF_LIBRARIES})
target_link_libraries(dw_bw ${TIFF_LIBRARIES})

#
# Gnu Scientific Library
#

target_link_libraries(dw gsl)
target_link_libraries(dw_bw gsl)

#
# OpenCL
#
# https://cmake.org/cmake/help/latest/module/FindOpenCL.html

if (ENABLE_GPU)
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
endif()


#
# OpenMP
#
#  https://cmake.org/cmake/help/v3.9/module/FindOpenMP.html
find_package(OpenMP)
if(OpenMP_C_FOUND)
  message("Found OpenMP")
  target_link_libraries(dw OpenMP::OpenMP_C)
  target_link_libraries(dw_bw OpenMP::OpenMP_C)
endif()


# Enable -flto (GCC) and similar if available, see
# https://cmake.org/cmake/help/latest/policy/CMP0069.html
check_ipo_supported(RESULT result)
if(result)
  message("Enabling IPO")
  set_property(DIRECTORY PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

include(GNUInstallDirs)
install(TARGETS dw)
install(TARGETS dw_bw)
install(FILES "${CMAKE_SOURCE_DIR}/doc/dw.1"
  DESTINATION "${CMAKE_INSTALL_MANDIR}/man1/" )
install(FILES "${CMAKE_SOURCE_DIR}/doc/dw_bw.1"
  DESTINATION "${CMAKE_INSTALL_MANDIR}/man1/" )
