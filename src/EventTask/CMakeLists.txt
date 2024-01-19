# Define the ParticleProperties library
add_library(EventTask_lib
        EventTask.cpp
        FreeMotionTask.cpp
)


# Specify the include directory for the ParticleProperties library
target_include_directories(EventTask_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
