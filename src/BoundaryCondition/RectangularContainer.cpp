#include "RectangularContainer.h"

// Implementation of RectangularContainer constructor
RectangularContainer::RectangularContainer(int id, const PropertyTypeID &type, int state,
                                           const Eigen::Vector3d &lowerCorner,
                                           const Eigen::Vector3d &dimensions,
                                           const Eigen::Vector3d &velocity)
    : BoundaryCondition(id, type, state), lowerCorner(lowerCorner), dimensions(dimensions)
{
    walls.resize(6);
    createWalls(velocity);
}

void RectangularContainer::rotateContainer(double angleDegrees, const Eigen::Vector3d &axis)
{
    // Convert angle to radians
    double angleRadians = angleDegrees * PI / 180.0;

    // Create a rotation matrix
    Eigen::AngleAxisd rotationMatrix(angleRadians, axis.normalized());

    // Rotate each wall
    for (auto &wall : walls)
    {
        // Rotate each corner point of the wall
        Eigen::Vector3d corner1Rotated = rotationMatrix * wall->getCorner1();
        Eigen::Vector3d corner2Rotated = rotationMatrix * wall->getCorner2();
        Eigen::Vector3d corner3Rotated = rotationMatrix * wall->getCorner3();
        Eigen::Vector3d normalRotated = rotationMatrix * wall->getNormal();

        // Update the wall with new corners (and possibly normal if needed)
        wall = std::make_shared<PlaneWall>(wall->getId(), wall->getType(), wall->getState(),
                         normalRotated, // Might need to rotate normal as well
                         corner1Rotated, corner2Rotated, corner3Rotated,
                         wall->getVelociy());
    }
}
void RectangularContainer::createWalls(const Eigen::Vector3d &velocity)
{
    // Assuming y is vertical, x and z are horizontal dimensions
    Eigen::Vector3d topNormal(0, -1, 0);   // Normal pointing upwards
    Eigen::Vector3d bottomNormal(0, 1, 0); // Normal pointing downwards
    Eigen::Vector3d leftNormal(1, 0, 0);   // Normal pointing left
    Eigen::Vector3d rightNormal(-1, 0, 0); // Normal pointing right
    Eigen::Vector3d frontNormal(0, 0, -1); // Normal pointing front
    Eigen::Vector3d backNormal(0, 0, 1);   // Normal pointing back

    // Create each wall with the appropriate corners and normals
    // You'll need to calculate the corners based on the lowerCorner and dimensions
    Eigen::Vector3d corner1 = lowerCorner;
    Eigen::Vector3d corner2 = lowerCorner + Eigen::Vector3d(0, 0, dimensions.z());
    Eigen::Vector3d corner3 = corner2 + Eigen::Vector3d(0, dimensions.y(), 0);
    Eigen::Vector3d corner4 = lowerCorner + Eigen::Vector3d(dimensions.x(), 0, 0);
    Eigen::Vector3d corner5 = corner4 + Eigen::Vector3d(0, 0, dimensions.z());
    Eigen::Vector3d corner6 = corner5 + Eigen::Vector3d(0, dimensions.y(), 0);
    Eigen::Vector3d corner7 = corner4 + Eigen::Vector3d(0, dimensions.y(), 0);
    Eigen::Vector3d corner8 = corner1 + Eigen::Vector3d(0, dimensions.y(), 0);

    walls[0] = std::make_shared<PlaneWall>(0, type, state, leftNormal, corner1, corner2, corner3, velocity);
    walls[1] = std::make_shared<PlaneWall>(1, type, state, rightNormal, corner4, corner5, corner6, velocity);
    walls[2] = std::make_shared<PlaneWall>(2, type, state, bottomNormal, corner4, corner5, corner2, velocity);
    walls[3] = std::make_shared<PlaneWall>(3, type, state, frontNormal, corner5, corner2, corner3, velocity);
    walls[4] = std::make_shared<PlaneWall>(4, type, state, backNormal, corner4, corner1, corner8, velocity);
    walls[5] = std::make_shared<PlaneWall>(5, type, state, topNormal, corner7, corner6, corner3, velocity);
}
std::string RectangularContainer::save_tostring() const
{
    std::ostringstream ss;
    ss << "RECTANGULARCONTAINER, " << id << ", " << type.getCategory() << ", " << type.getSubType() << ", " << state << ", "
       << lowerCorner.x() << ", " << lowerCorner.y() << ", " << lowerCorner.z() << ", "
       << dimensions.x() << ", " << dimensions.y() << ", " << dimensions.z() << ", "
       << walls[0]->getVelociy().x() << ", " << walls[0]->getVelociy().y() << ", " << walls[0]->getVelociy().z();
    return ss.str();
}