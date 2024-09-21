target_sources(geodesy 
  PRIVATE
    ${CMAKE_SOURCE_DIR}/src/cartesian_to_spherical.cpp  
    ${CMAKE_SOURCE_DIR}/src/geodetic_to_lvlh.cpp
    ${CMAKE_SOURCE_DIR}/src/spherical_to_cartesian.cpp
)
