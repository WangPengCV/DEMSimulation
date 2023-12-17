cmake_minimum_required(VERSION 3.10)
project(DEMSimulation)


set(EIGEN3_INCLUDE_DIR "D:/opensource/Eigen/Eigen_install/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})

# Add src and tests subdirectories
add_subdirectory(src)
add_subdirectory(tests)