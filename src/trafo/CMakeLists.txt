cmake_minimum_required(VERSION 3.22)

project(trafo
  LANGUAGES C)

set (CMAKE_C_STANDARD 11)

if(WIN32)
  add_definitions(-DWINDOWS)
endif()

add_definitions(-D_DEFAULT_SOURCE)
option(BUILD_SHARED_LIBS "Build shared library" ON)
# To build static:  -DBUILD_SHARED_LIBS=OFF

set(CMAKE_C_VISIBILITY_PRESET "hidden")
set(CMAKE_VISIBILITY_INLINES_HIDDEN TRUE)

#
# Default build type is RELEASE
#
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if (NOT CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DNDEBUG)
endif()

add_library(trafo
  src/entropy.c
  src/ftab.c
  src/gini.c
  src/sortbox.c
  src/trafo.c
  src/trafo_util.c
  src/qsort.c)

add_executable(trafo_cli
  src/trafo_cli.c
  src/entropy.c
  src/ftab.c
  src/gini.c
  src/qsort.c
  src/sortbox.c
  src/trafo.c
  src/trafo_cli.c
  src/trafo_util.c)

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
INSTALL(TARGETS trafo
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
