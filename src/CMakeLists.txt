# Include the ParticleProperties directory to build the ParticleProperties_lib
add_subdirectory(ParticleProperties)
add_subdirectory(Particle)
add_subdirectory(BoundaryCondition)
add_subdirectory(Visualization)

# Define the DEM_lib library, which may include other source files in the src directory
add_library(DEM_lib
    DEMProperties.cpp
    DEMProperties.h
)

# Link the DEM_lib with ParticleProperties_lib if needed
target_link_libraries(DEM_lib
    ParticleProperties_lib
    Particle_lib
    BoundaryCondition_lib
    Visualization_lib
)

# Specify the include directories for the DEM_lib library
target_include_directories(DEM_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/ParticleProperties
    ${CMAKE_CURRENT_SOURCE_DIR}/Particle
    ${CMAKE_CURRENT_SOURCE_DIR}/BoundaryCondition
    ${CMAKE_CURRENT_SOURCE_DIR}/Visualization
)
