# Define the ParticleProperties library
add_library(ContactDetection_lib
    GridBasedContactDetection.cpp
)

# Link with ParticleProperties_lib to use PropertyTypeID
target_link_libraries(ContactDetection_lib
    Particle_lib
    ParticleProperties_lib
    BoundaryCondition_lib
)
# Specify the include directory for the ParticleProperties library
target_include_directories(ContactDetection_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
