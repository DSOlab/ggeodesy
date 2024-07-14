cmake_minimum_required(VERSION 3.1)

project(
  Geodesy
  VERSION 1.1.0
  DESCRIPTION "DSO Geodetic Library"
  LANGUAGES CXX)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  
  # Let's ensure -std=c++xx instead of -std=g++xx
  set(CMAKE_CXX_EXTENSIONS OFF)
  
  # Let's nicely support folders in IDEs
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
  
  include(CTest)
  
  find_package(Doxygen)
  if(Doxygen_FOUND)
    add_subdirectory(doc)
  else()
    message(STATUS "Doxygen not found, not building docs")
  endif()

endif()

# We need Eigen
find_package(Eigen3 3.3 REQUIRED)

# The library
add_subdirectory(src)

# The tests
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
