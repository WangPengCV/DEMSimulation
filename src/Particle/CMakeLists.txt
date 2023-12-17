# Define the ParticleProperties library
add_library(Particle_lib
    Particle.cpp
    SphereParticle.cpp
)

# Link with ParticleProperties_lib to use PropertyTypeID
target_link_libraries(Particle_lib
    ParticleProperties_lib
)
# Specify the include directory for the ParticleProperties library
target_include_directories(Particle_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
