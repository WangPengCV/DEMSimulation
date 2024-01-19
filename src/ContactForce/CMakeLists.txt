# Define the ParticleProperties library
add_library(ContactForce_lib
    ContactForce.cpp
)

# Link with ParticleProperties_lib to use PropertyTypeID
target_link_libraries(ContactForce_lib
    Particle_lib
    ParticleProperties_lib
    BoundaryCondition_lib
)
# Specify the include directory for the ParticleProperties library
target_include_directories(ContactForce_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
