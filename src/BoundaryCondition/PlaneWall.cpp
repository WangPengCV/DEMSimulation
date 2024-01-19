#include "PlaneWall.h"

<<<<<<< HEAD
=======
#include "PlaneWall.h"

>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
void PlaneWall::generateMesh(double meshResolution)
{
    meshVertices.clear(); // Clear existing vertices

<<<<<<< HEAD
    // Calculate the basis vectors for the plane
=======
   // Calculate the basis vectors for the plane
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    Eigen::Vector3d u = (corner2 - corner1).normalized();
    Eigen::Vector3d v = (corner3 - corner2).normalized();

    // Determine the number of divisions along each basis vector
    double lengthU = (corner2 - corner1).norm();
    double lengthV = (corner3 - corner2).norm();
    int divisionsU = static_cast<int>(lengthU / meshResolution);
    int divisionsV = static_cast<int>(lengthV / meshResolution);

<<<<<<< HEAD
    // Generate grid points
    for (int i = 0; i <= divisionsU; ++i)
    {
        for (int j = 0; j <= divisionsV; ++j)
        {
=======
   
    // Generate grid points
    for (int i = 0; i <= divisionsU; ++i) 
    {
        for (int j = 0; j <= divisionsV; ++j) {
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
            Eigen::Vector3d point = corner1 + u * (meshResolution * i) + v * (meshResolution * j);
            meshVertices.push_back(point);
        }
    }
}
<<<<<<< HEAD
=======

>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
const std::vector<Eigen::Vector3d> &PlaneWall::getMeshVertices() const
{
    return meshVertices;
}
<<<<<<< HEAD
=======

>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
const Eigen::Vector3d &PlaneWall::getNormal() const
{
    return normal;
}
const Eigen::Vector3d &PlaneWall::getCorner1() const
{
    return corner1;
}
const Eigen::Vector3d &PlaneWall::getCorner2() const
{
    return corner2;
}
const Eigen::Vector3d &PlaneWall::getCorner3() const
{
    return corner3;
}
<<<<<<< HEAD
const Eigen::Vector3d &PlaneWall::getVelocity() const
{
    return velocity;
}
void PlaneWall::setNormal(Eigen::Vector3d &Normal)
{
    normal = Normal;
}
void PlaneWall::setCorner1(Eigen::Vector3d &Corner1)
{
    corner1 = Corner1;
}
void PlaneWall::setCorner2(Eigen::Vector3d &Corner2)
{
    corner2 = Corner2;
}
void PlaneWall::setCorner3(Eigen::Vector3d &Corner3)
{
    corner3 = Corner3;
}
void PlaneWall::setVelocity(Eigen::Vector3d &Velociy)
{
    velocity = Velociy;
=======

const Eigen::Vector3d &PlaneWall::getVelociy() const
{
    return velocity;
}

void PlaneWall::setNormal(Eigen::Vector3d &normal)
{
    normal = normal;
}
void PlaneWall::setCorner1(Eigen::Vector3d &corner1)
{
    corner1 = corner1;
}
void PlaneWall::setCorner2(Eigen::Vector3d &corner2)
{
    corner2 = corner2;
}
void PlaneWall::setCorner3(Eigen::Vector3d &corner3)
{
    corner3 = corner3;
}
void PlaneWall::setVelociy(Eigen::Vector3d &velociy)
{
    velociy = velociy;
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
}
void PlaneWall::addForce(const Eigen::Vector3d &additionalForce)
{
    force += additionalForce;
<<<<<<< HEAD
}
void PlaneWall::resetForce()
{
    force.setZero();
}

std::string PlaneWall::save_tostring() const
{
=======
}   
void PlaneWall::resetForce()
{
        force.setZero();

}

std::string PlaneWall::save_tostring() const {
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    std::ostringstream ss;
    ss << "PLANEWALL, " << id << ", " << type.getCategory() << ", " << type.getSubType() << ", " << state << ", "
       << normal.x() << ", " << normal.y() << ", " << normal.z() << ", "
       << corner1.x() << ", " << corner1.y() << ", " << corner1.z() << ", "
       << corner2.x() << ", " << corner2.y() << ", " << corner2.z() << ", "
       << corner3.x() << ", " << corner3.y() << ", " << corner3.z() << ", "
       << velocity.x() << ", " << velocity.y() << ", " << velocity.z();
    return ss.str();
}