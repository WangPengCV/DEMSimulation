#pragma once
#include "BoundaryCondition.h"

class PlaneWall : public BoundaryCondition
{
public:
    PlaneWall(int id, const PropertyTypeID &type, int state, const Eigen::Vector3d &normal,
              const Eigen::Vector3d &corner1, Eigen::Vector3d& corner2,Eigen::Vector3d& corner3,double meshResolution)
        : BoundaryCondition(id, type, state), normal(normal), corner1(corner1),corner2(corner2),corner3(corner3),meshResolution(meshResolution) {}

    const std::vector<Eigen::Vector3d> &getMeshVertices() const;
    const Eigen::Vector3d &getNormal() const;
    const Eigen::Vector3d &getCorner1() const;
    const Eigen::Vector3d &getCorner2() const;
    const Eigen::Vector3d &getCorner3() const;


private:
    void generateMesh();

    Eigen::Vector3d normal;     // Normal to the plane
    Eigen::Vector3d corner1, corner2, corner3; //Three corner points of a plane wall (arranged clockwise)  
    std::vector<Eigen::Vector3d> meshVertices; // discreted points for contact detection
    double meshResolution; // mesh size smaller than grid size

};
