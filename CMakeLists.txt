# Usage:
# mkdir build
# cd build
# cmake ..
# or
# CC=clang cmake ..
# make dw

# See https://medium.com/@onur.dundar1/cmake-tutorial-585dd180109b
# https://cliutils.gitlab.io/modern-cmake/

cmake_minimum_required(VERSION 3.9)
project(deconwolf VERSION 0.3.6
  DESCRIPTION "Deconvolution of wide field images"
  LANGUAGES C)

if(${CMAKE_SYSTEM_NAME} MATCHES Darwin)
  message("This is a ${CMAKE_SYSTEM_NAME} system")
elseif(${CMAKE_SYSTEM_NAME} MATCHES Windows)
  message("This is a Windows System")
endif()

# UNIX, WIN32, WINRT, CYGWIN, APPLE are environment variables as flags set by default system
if(UNIX)
  message("This is a ${CMAKE_SYSTEM_NAME} system")
elseif(WIN32)
  message("This is a Windows System")
endif()

add_definitions(-DMyDefine="Windows PC")





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
  src/method_eve.c
  src/method_identity.c
  src/method_rl.c
  src/method_ave.c
  src/method_eve.c
  src/method_shb.c
  src/tiling.c
  src/fft.c)


target_link_libraries(dw fftw3 fftw3f fftw3f_omp)

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
  target_link_libraries(dw ${MATH_LIBRARY})
endif()

target_link_libraries(dw tiff)
target_link_libraries(dw gsl)

# https://cmake.org/cmake/help/v3.9/module/FindOpenMP.html
find_package(OpenMP)
if(OpenMP_C_FOUND)
  message("Found OpenMP")
  message(${OpenMP_C_FLAGS}) # how to set these flags later on?
  message(${OpenMP_C_LIB_NAMES}) #

  target_link_libraries(dw ${OpenMP_C_LIB_NAMES})
endif()
