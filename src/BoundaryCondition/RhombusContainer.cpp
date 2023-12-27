#include "RhombusContainer.h"

RhombusContainer::RhombusContainer(int id, const PropertyTypeID& type, int state, 
                                   const Eigen::Vector3d& corner1, const Eigen::Vector3d& corner2,
                                   const Eigen::Vector3d& corner3, const Eigen::Vector3d& corner4,
                                   const Eigen::Vector3d& velocity)
    : BoundaryCondition(id, type, state) {
    // Set up the walls based on the corners
    createWalls(velocity);
}

void RhombusContainer::createWalls(const Eigen::Vector3d& velocity) {
    // Calculate the normal vectors for each wall and create PlaneWall objects
    // Assuming the corners are provided in a consistent order
    // Example for one wall:
    Eigen::Vector3d normal = (corner2 - corner1).cross(corner3 - corner1).normalized();
    walls[0] = PlaneWall(/* parameters for wall 0 */, normal, corner1, corner2, corner3, velocity);
    // Repeat for other walls
}

std::string RhombusContainer::save_tostring() const {
    // Implement string representation for saving
}

PlaneWall& RhombusContainer::getWall(int index) {
    return walls[index];
}
