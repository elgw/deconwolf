cmake_minimum_required(VERSION 3.9)

project(kdtree
  LANGUAGES C)

set (CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(BUILD_SHARED_LIBS "Build shared library" ON)
# To build static:  -DBUILD_SHARED_LIBS=OFF
set(CMAKE_C_VISIBILITY_PRESET "hidden")
set(CMAKE_VISIBILITY_INLINES_HIDDEN TRUE)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


add_library( kdtree src/kdtree.c src/pqheap.c src/quickselect.c )



target_include_directories(kdtree PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

set_target_properties(kdtree PROPERTIES PUBLIC_HEADER "include/kdtree.h")

include(GNUInstallDirs)

INSTALL(TARGETS kdtree
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
