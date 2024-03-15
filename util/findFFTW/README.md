CMake module for finding FFTW 3 using find_package

# Usage

Once added to your project, this module allows you to find FFTW libraries and headers using the CMake `find_package` command:

```cmake
find_package(FFTW [REQUIRED] [QUIET] [COMPONENTS component1 ... componentX] )
```

This module sets the following variables:
- `FFTW_FOUND`                  ... true if fftw is found on the system
- `FFTW_[component]_LIB_FOUND`  ... true if the component is found on the system (see components below)
- `FFTW_LIBRARIES`              ... full paths to all found fftw libraries
- `FFTW_[component]_LIB`        ... full path to one of the components (see below)
- `FFTW_INCLUDE_DIRS`           ... fftw include directory paths

The following variables will be checked by the module:
- `FFTW_USE_STATIC_LIBS`        ... if true, only static libraries are found, otherwise both static and shared.
- `FFTW_ROOT`                   ... if set, the libraries are exclusively searched under this path.

This package supports the following components:
- `FLOAT_LIB`
- `DOUBLE_LIB`
- `LONGDOUBLE_LIB`
- `FLOAT_THREADS_LIB`
- `DOUBLE_THREADS_LIB`
- `LONGDOUBLE_THREADS_LIB`
- `FLOAT_OPENMP_LIB`
- `DOUBLE_OPENMP_LIB`
- `LONGDOUBLE_OPENMP_LIB`
- `FLOAT_MPI_LIB`
- `DOUBLE_MPI_LIB`
- `LONGDOUBLE_MPI_LIB`

and the following linking targets

- `FFTW::Float`
- `FFTW::Double`
- `FFTW::LongDouble`
- `FFTW::FloatThreads`
- `FFTW::DoubleThreads`
- `FFTW::LongDoubleThreads`
- `FFTW::FloatOpenMP`
- `FFTW::DoubleOpenMP`
- `FFTW::LongDoubleOpenMP`
- `FFTW::FloatMPI`
- `FFTW::DoubleMPI`
- `FFTW::LongDoubleMPI`

# Adding to your project

## Automatic download from CMake project

Copy the following into the `CMakeLists.txt` file of the project you want to use FindFFTW in:
```cmake
configure_file(downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
    else()
    message("CMake step for findFFTW completed (${result}).")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif()

set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")
```

And add a file called `downloadFindFFTW.cmake.in` to your project containing the following:
```cmake
cmake_minimum_required(VERSION 2.8.2)

project(findFFTW-download NONE)

include(ExternalProject)

ExternalProject_Add(findFFTW_download
    GIT_REPOSITORY    "https://github.com/egpbos/findfftw.git"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    TEST_COMMAND      ""
    SOURCE_DIR        "${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src"
    BINARY_DIR        ""
    INSTALL_DIR       ""
)
```

After this, `find_package(FFTW)` can be used in the `CMakeLists.txt` file.

## Manual

Clone the repository into directory `PREFIX/findFFTW`:
```sh
git clone https://github.com/egpbos/findfftw.git PREFIX/findFFTW
```

Then add the following to your `CMakeLists.txt` to allow CMake to find the module:
```cmake
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "PREFIX/findFFTW")
```
