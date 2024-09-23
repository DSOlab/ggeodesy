cmake_minimum_required(VERSION 3.11)

project(
  Geodesy
  VERSION 1.1.0
  DESCRIPTION "DSO Geodetic Library"
  LANGUAGES CXX)

# We need Eigen
find_package(Eigen3 3.3 REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)
add_compile_options(-Wall 
  -Wextra 
  -Werror 
  -pedantic 
  -W 
  -Wshadow 
  $<$<CONFIG:Release>:-O2>
  $<$<CONFIG:Release>:-march=native>
  $<$<CONFIG:Debug>:-g>
  $<$<CONFIG:Debug>:-pg> 
  $<$<CONFIG:Debug>:-Wdisabled-optimization>
  )
add_compile_definitions(
  $<$<CONFIG:Debug>:DEBUG>
)

add_library(geodesy)
target_include_directories(geodesy PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
  $<INSTALL_INTERFACE:include/geodesy>
)

add_subdirectory(src)

# The tests
include(CTest)
add_subdirectory(test/unit_tests)
enable_testing()

# Install headers at: $PREFIX/datetime/...
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
	DESTINATION include/geodesy
)

# install library
install(TARGETS geodesy
        EXPORT geodesyTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
)

install(EXPORT geodesyTargets
        FILE geodesyTargets.cmake
        NAMESPACE dso::
        DESTINATION lib/cmake/geodesy
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "geodesyConfigVersion.cmake"
  VERSION ${geodesy_version}
  COMPATIBILITY AnyNewerVersion
)

install(FILES "geodesyConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/geodesyConfigVersion.cmake"
        DESTINATION lib/cmake/dso
)
