cmake_minimum_required(VERSION 3.10)
project(DEMSimulation)




set(VTK_DIR "D:/opensource/vtk/build")
find_package(VTK REQUIRED)


set(EIGEN3_INCLUDE_DIR "D:/opensource/Eigen/Eigen_install/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(${PROJECT_NAME} main.cpp Simulation.cpp Visualization.cpp Contactforce.cpp Spatialgrid.cpp Datatype.h Particle.h)
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} )

vtk_module_autoinit(
    TARGETS ${PROJECT_NAME}
    MODULES ${VTK_LIBRARIES}
  )