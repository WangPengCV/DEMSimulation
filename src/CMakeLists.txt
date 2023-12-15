# Create a library from source files
add_library(particle_lib ParticleProperties.cpp ParticlePropertyManager.cpp DEMProperties.cpp)

# Make sure the compiler can find include files for our library
# when other libraries or executables link to particle_lib
target_include_directories(particle_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
