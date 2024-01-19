# Define the ParticleProperties library
add_library(Particle_lib
    Particle.cpp
    SphereParticle.cpp
<<<<<<< HEAD
    SphereCylinderBond.cpp
    Fiber.cpp
=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
)

# Link with ParticleProperties_lib to use PropertyTypeID
target_link_libraries(Particle_lib
    ParticleProperties_lib
    BoundaryCondition_lib
    ContactForce_lib
)
# Specify the include directory for the ParticleProperties library
target_include_directories(Particle_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
