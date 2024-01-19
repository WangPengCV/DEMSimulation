# Define the ParticleProperties library
add_library(ParticleProperties_lib
    ParticleProperties.cpp
    SphereProperties.cpp
    PlanewallProperties.cpp
    ParticlePropertyManager.cpp
<<<<<<< HEAD
    FiberProperties.cpp
=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
)

# Specify the include directory for the ParticleProperties library
target_include_directories(ParticleProperties_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
