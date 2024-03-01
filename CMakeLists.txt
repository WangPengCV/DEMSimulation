cmake_minimum_required(VERSION 3.12)
# Set C++17 as the standard for the whole project
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(DEMSimulation)

set(VTK_DIR "D:/opensource/vtk/build")
find_package(VTK REQUIRED)
include_directories(${VTK_INCLUDE_DIRS})

set(EIGEN3_INCLUDE_DIR "D:/opensource/Eigen/Eigen_install/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})

set(MPI_INCLUDE_PATH "D:/opensource/MPI/install/Include")
set(MPI_CXX_LIBRARIES "D:/opensource/MPI/install/Lib/x64/msmpi.lib") # 假设是64位库
include_directories(${MPI_INCLUDE_PATH})


# Add src and tests subdirectories
add_subdirectory(src)
add_subdirectory(tests)

