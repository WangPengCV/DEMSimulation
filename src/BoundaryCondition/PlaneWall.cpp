#include "PlaneWall.h"

#include "PlaneWall.h"

void PlaneWall::generateMesh()
{
    meshVertices.clear(); // Clear existing vertices

   // Calculate the basis vectors for the plane
    Eigen::Vector3d u = (corner2 - corner1).normalized();
    Eigen::Vector3d v = (corner3 - corner2).normalized();

    // Determine the number of divisions along each basis vector
    double lengthU = (corner2 - corner1).norm();
    double lengthV = (corner3 - corner2).norm();
    int divisionsU = static_cast<int>(lengthU / meshResolution);
    int divisionsV = static_cast<int>(lengthV / meshResolution);

   
    // Generate grid points
    for (int i = 0; i <= divisionsU; ++i) 
    {
        for (int j = 0; j <= divisionsV; ++j) {
            Eigen::Vector3d point = corner1 + u * (meshResolution * i) + v * (meshResolution * j);
            meshVertices.push_back(point);
        }
    }
}

const std::vector<Eigen::Vector3d> &PlaneWall::getMeshVertices() const
{
    return meshVertices;
}

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
