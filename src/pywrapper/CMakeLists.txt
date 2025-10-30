cmake_minimum_required(VERSION 3.28 FATAL_ERROR)
project(deconwolf_python
  VERSION 1.0
  DESCRIPTION "Python module for deconwolf"
  LANGUAGES C)

set (CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Happy to load libdeconwolf from the same folder.
set(CMAKE_INSTALL_RPATH "./")

# Prefer that Python use whatever environment that
# is active. Especially the system version of Numpy can old.
set (Python_FIND_VIRTUALENV FIRST) # Note: Python, not Python3 here
find_package(Python3 3.8 REQUIRED
  COMPONENTS Development.Module NumPy)

Python3_add_library(pydeconwolf src/pydeconwolf.c)
# link do libdeconwolf[.so/.dll] which needs be in place to use the module
target_link_libraries(pydeconwolf PRIVATE deconwolf Python3::NumPy)
