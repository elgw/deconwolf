cmake_minimum_required(VERSION 3.22)

project(trafo
  LANGUAGES C)

set (CMAKE_C_STANDARD 11)

if(WIN32)
  add_definitions(-DWINDOWS)
endif()

#
# Default build type is RELEASE
#
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Don't install to /usr/local/ use /usr/
set( CMAKE_INSTALL_PREFIX
  /usr )

add_library(trafo SHARED
  src/entropy.c
  src/ftab.c
  src/gini.c
  src/sortbox.c
  src/trafo.c
  src/trafo_util.c
  src/qsort.c)


set(TRAFO_PUBLIC_HEADERS
  "src/trafo.h"
)

add_executable(trafo_cli
  src/entropy.c
  src/ftab.c
  src/gini.c
  src/qsort.c
  src/sortbox.c
  src/trafo.c
  src/trafo_cli.c
  src/trafo_util.c)

target_link_libraries(trafo_cli
  trafo)

set_target_properties(trafo_cli
  PROPERTIES INSTALL_RPATH "$ORIGIN")

set_target_properties(trafo
  PROPERTIES PUBLIC_HEADER src/trafo.h)

find_library(MATH_LIBRARY m)
if(MATH_LIBRARY)
  target_link_libraries(trafo ${MATH_LIBRARY})
  target_link_libraries(trafo_cli ${MATH_LIBRARY})
endif()

#
# Getopt
#
if(WIN32)
  find_library(GETOPT_LIBRARY getopt)
  message("getopt found at: ${GETOPT_LIBRARY}")
  target_link_libraries(trafo_cli ${GETOPT_LIBRARY})
  SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GETOPT_LIBRARY}")
  find_path(GETOPT_INCLUDE_DIR getopt.h)
  target_include_directories(trafo_cli PRIVATE GETOPT_INCLUDE_DIR)
endif()


#
# OpenMP
#
#  https://cmake.org/cmake/help/v3.9/module/FindOpenMP.html
find_package(OpenMP)
if(OpenMP_C_FOUND)
  message("Found OpenMP")
  target_link_libraries(trafo OpenMP::OpenMP_C)
  target_link_libraries(trafo_cli OpenMP::OpenMP_C)
endif()

#
# Only export symbols explicitly marked for export.
#
# This should work with gcc/linux and could probably be adjusted for
# other compiler/os as well.
#

include(CheckCSourceCompiles)
check_c_source_compiles(
  "
        #include <stdlib.h>
        static void f(void) __attribute__ ((visibility(\"hidden\")));
        int main(void) {return 0;}
    "
  HAVE___ATTRIBUTE__VISIBILITY_HIDDEN
)

if(HAVE___ATTRIBUTE__VISIBILITY_HIDDEN)
  # target_compile_options(trafo PRIVATE -DHAVE___ATTRIBUTE__VISIBILITY_HIDDEN -fvisibility=hidden)
  target_compile_options(trafo PRIVATE -DHAVE___ATTRIBUTE__VISIBILITY_HIDDEN)
  set_target_properties(trafo PROPERTIES C_VISIBILITY_PRESET hidden)
endif()

#
# Link time optimization
#

# https://cmake.org/cmake/help/latest/module/CheckIPOSupported.html
include(CheckIPOSupported)
check_ipo_supported(RESULT result OUTPUT output)
if(result)
  message("Enabling IPO")
  set_property(TARGET trafo PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set_property(TARGET trafo_cli PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_JOBS 8)
else()
  message(WARNING "IPO is not supported: ${output}")
endif()

include("GNUInstallDirs")

set_target_properties( trafo PROPERTIES
  PUBLIC_HEADER "${TRAFO_PUBLIC_HEADERS}")

INSTALL(TARGETS trafo
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_FULL_INCLUDEDIR}"
  LIBRARY DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
)
