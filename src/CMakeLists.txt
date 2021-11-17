cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(deconwolf LANGUAGES C)

set(WINDOWS, FALSE)
#set(WINDOWS TRUE)

set(CMAKE_C_FLAGS_DEBUG_INIT "-Wall -O0 -g -std=gnu99")
set(CMAKE_C_FLAGS_RELEASE_INIT "-Wall -O3 -DNDEBUG -flto -march=native -std=gnu99")

INCLUDE_DIRECTORIES(/usr/include)
LINK_DIRECTORIES(/usr/lib/x86_64-linux-gnu)

add_executable(dw
  dw.c
  fft.c
  fim.c
  fim_tiff.c
  deconwolf.c
  tiling.c
  deconwolf_tif_max.c
  )

add_executable(dwbw
  dw_bwpsf.c
  lanczos3.c
  fim.c
  fft.c
  fim_tiff.c
  bw_gsl.c)

if(${WINDOWS})
  set(CMAKE_C_FLAGS "-DWINDOWS -Wall -Wextra -g")
  set(TIFF_LIBRARY tiff)
  set(TIFF_INCLUDE_DIR ../win/libtiff/)
  set(TIFF_LIB_DIR ../win/)
  set(FFTW3f_LIBRARY fftw3f-3)
  set(FFTW3f_INCLUDE_DIR ../win/fftw3/)
  set(GSL_LIBRARY libgsl-25)
  set(GSL_INCLUDE_DIR ../win/)
else()
  set(TIFF_LIBRARY /usr/lib/x86_64-linux-gnu/libtiff.so)
  set(TIFF_INCLUDE_DIR /usr/include/x86_64-linux-gnu)
  set(FFTW3f_LIBRARY /usr/lib/x86_64-linux-gnu/libfftw3f.so)
  set(FFTW3f_INCLUDE_DIR /usr/include)
  set(FFTW3ft_LIBRARY /usr/lib/x86_64-linux-gnu/libfftw3f_threads.so)
  set(FFTW3ft_INCLUDE_DIR /usr/include)
endif()

find_package(GSL REQUIRED)

include(CheckLibraryExists)
CHECK_LIBRARY_EXISTS(m sin "" HAVE_LIB_M)
CHECK_LIBRARY_EXISTS($(FFTW3f_LIBRARY) fftwf_execute "" HAVE_LIB_FFTW3)

target_link_libraries(dw
  ${TIFF_LIBRARY}
  ${FFTW3f_LIBRARY}
  ${FFTW3ft_LIBRARY}
  ${GSL_LIBRARY}
  m)


target_include_directories(dw PRIVATE
  ${TIFF_INCLUDE_DIR}
  ${FFTW3f_INCLUDE_DIR}
  ${GSL_INCLUDE_DIR})

target_link_libraries(dwbw
  ${TIFF_LIBRARY}
  ${GSL_LIBRARY}
  ${FFTW3f_LIBRARY}
  ${FFTW3ft_LIBRARY}
  ${GSL_LIBRARY}
  m)

target_include_directories(dwbw PRIVATE
  ${TIFF_INCLUDE_DIR}
  ${FFTW3f_INCLUDE_DIR}
  ${GSL_INCLUDE_DIR}
  gsl)

target_link_directories(dwbw
  PUBLIC ${TIFF_LIB_DIR})

target_link_directories(dw
  PUBLIC ${TIFF_LIB_DIR})

if(${WINDOWS})
  #set(CMAKE_C_COMPILER i686-w64-mingw32-gcc)
  set(CMAKE_C_COMPILER x86_64-w64-mingw32-gcc)
else()
  set(CMAKE_C_COMPILER gcc)
endif()
