#pragma once
#include "BoundaryCondition.h"
#include "PlaneWall.h"
#include <array>
#include <Eigen/Dense>

class RhombusContainer : public BoundaryCondition {
public:
    // Constructor
    RhombusContainer(int id, const PropertyTypeID& type, int state, 
                     const Eigen::Vector3d& corner1, const Eigen::Vector3d& corner2,
                     const Eigen::Vector3d& corner3, const Eigen::Vector3d& corner4,
                     const Eigen::Vector3d& velocity = Eigen::Vector3d::Zero());

    // Override the save_tostring method from BoundaryCondition
    virtual std::string save_tostring() const override;

    // Accessors for walls
    PlaneWall& getWall(int index);

private:
    std::array<PlaneWall, 4> walls;

    // Helper method to create walls
    void createWalls(const Eigen::Vector3d& velocity);
};
