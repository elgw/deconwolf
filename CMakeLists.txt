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
#   Either install directly with
#   cmake install
#
#   Or produce a package that can be installed with the
#   system package manager, for example
#   cpack -G RPM
#   cpack -G DEB
#
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
# cmake documentation and tutorials:
# https://cmake.org/cmake/help/latest/
# https://cliutils.gitlab.io/modern-cmake/
#

# TODO:
# Check https://learn.microsoft.com/en-us/cpp/error-messages/compiler-warnings/compiler-warning-level-3-c4996?view=msvc-170

cmake_minimum_required(VERSION 3.9)

project(deconwolf
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
  src/method_identity.c
  src/method_rl.c
  src/method_shb.c
  src/tiling.c
  src/fft.c
  src/gmlfit.c
  src/sparse_preprocess.c
  src/sparse_preprocess_cli.c)

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
# Parse version number from src/dw_version.h
#

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/src/dw_version.h DW_VERSION_MAJOR
  REGEX "#define[ ]+DW_VERSION_MAJOR[ ]+\"[0-9]+\"")
string(REGEX MATCH "\"([0-9]+)\"" _ ${DW_VERSION_MAJOR})
set(DW_VERSION_MAJOR "${CMAKE_MATCH_1}")

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/src/dw_version.h DW_VERSION_MINOR
  REGEX "#define[ ]+DW_VERSION_MINOR[ ]+\"[0-9]+\"")
string(REGEX MATCH "\"([0-9]+)\"" _ ${DW_VERSION_MINOR})
set(DW_VERSION_MINOR "${CMAKE_MATCH_1}")

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/src/dw_version.h DW_VERSION_PATCH
  REGEX "#define[ ]+DW_VERSION_PATCH[ ]+\"[0-9]+\"")
string(REGEX MATCH "\"([0-9]+)\"" _ ${DW_VERSION_PATCH})
set(DW_VERSION_PATCH "${CMAKE_MATCH_1}")

project(deconwolf VERSION "${DW_VERSION_MAJOR}.${DW_VERSION_MINOR}.${DW_VERSION_PATCH}")

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
# Getopt
#
if(WIN32)
# find_path(GETOPT_INCLUDE_DIR getopt.h)
find_library(GETOPT_LIBRARY getopt)
message("getopt found at: ${GETOPT_LIBRARY}")
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GETOPT_LIBRARY}")
endif()

#
# FFTW3
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/util/findFFTW)
if(UNIX)
find_package(FFTW REQUIRED FLOAT_LIB FLOAT_OPENMP_LIB)
else()
find_package(FFTW REQUIRED FLOAT_LIB)
endif()
target_link_libraries(dw FFTW::Double)
target_link_libraries(dw FFTW::Float)


target_link_libraries(dw_bw FFTW::Float)

if(UNIX)
target_link_libraries(dw_bw FFTW::FloatOpenMP)
target_link_libraries(dw FFTW::FloatOpenMP)
endif()


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
find_package(GSL)
target_link_libraries(dw ${GSL_LIBRARIES})
target_link_libraries(dw_bw ${GSL_LIBRARIES})

#
# OpenCL
#
# https://cmake.org/cmake/help/latest/module/FindOpenCL.html

if (ENABLE_GPU)
find_package(OpenCL)
if(OpenCL_FOUND)
  message("Will use OpenCL")
  target_link_libraries(dw ${OpenCL_LIBRARIES})
  target_include_directories(dw PRIVATE ${OpenCL_INCLUDE_DIRS})
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
else()
message("OpenMP not found automatically, trying manual")

set(
  CMAKE_C_FLAGS
  "${CMAKE_C_FLAGS} -fopenmp=omp /openmp"
)
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fopenmp=omp -lomp")

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

# https://cmake.org/cmake/help/book/mastering-cmake/chapter/Packaging%20With%20CPack.html


SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Erik Wernersson")
# resolve dependencies automatically
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "https://www.github.com/elgw/deconwolf/")
SET(CPACK_DEBIAN_PACKAGE_SECTION "science")
INCLUDE(CPack)

# The produced deb files are not perfectly portable,
# dpkg -I ./deconwolf-0.3.8-Linux.deb
# reveals ocl-icd-libopencl1-amdgpu-pro
