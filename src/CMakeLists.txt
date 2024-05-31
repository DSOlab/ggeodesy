add_library(geodesy)

target_sources(geodesy 
  PRIVATE
    cartesian_to_spherical.cpp  
    geodetic_to_lvlh.cpp
    spherical_to_cartesian.cpp
)

target_include_directories(geodesy
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../include/core
)

# We need the Eigen-3 library
target_link_libraries(geodesy PRIVATE Eigen3::Eigen)

# Install as: $PREFIX/geodesy/...
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include/
	DESTINATION geodesy
)
