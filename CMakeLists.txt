cmake_minimum_required(VERSION 3.1)

project(
  Geodesy
  VERSION 1.1.0
  DESCRIPTION "DSO Geodetic Library"
  LANGUAGES CXX)

# We need Eigen
find_package(Eigen3 3.3 REQUIRED)

# The library
add_subdirectory(src)

# The tests
include(CTest)
add_subdirectory(test//unit_tests)
enable_testing()

add_compile_options(
  -Wall -Wextra -Werror -pedantic -W -Wshadow -march=native -DEIGEN_NO_AUTOMATIC_RESIZING
  $<$<CONFIG:RELEASE>:-O2 -march=native -DEIGEN_NO_AUTOMATIC_RESIZING>
  $<$<CONFIG:DEBUG>:-g -pg -Wdisabled-optimization -DDEBUG>
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)
