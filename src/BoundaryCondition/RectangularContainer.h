#pragma once
#include "BoundaryCondition.h"
#include "PlaneWall.h"
#include <unordered_map>
#include <Eigen/Dense>

class RectangularContainer : public BoundaryCondition {
public:
    RectangularContainer();
    // Constructor
    RectangularContainer(int id, const PropertyTypeID& type, int state, 
                         const Eigen::Vector3d& lowerCorner, 
                         const Eigen::Vector3d& dimensions, 
                         const Eigen::Vector3d& velocity = Eigen::Vector3d::Zero());

    void rotateContainer(double angleDegrees, const Eigen::Vector3d& axis);

    const std::vector<std::shared_ptr<PlaneWall>>& getPlaneWall() const {return walls;} 

    // Override the save_tostring method from the BoundaryCondition class
    virtual std::string save_tostring() const override;

    // Methods to access individual walls
    const std::shared_ptr<PlaneWall>& getTopWall() const { return walls[5];}
    const std::shared_ptr<PlaneWall>& getBottomWall() const {return walls[2];} 
    const std::shared_ptr<PlaneWall>& getLeftWall() const {return walls[0];}
    const std::shared_ptr<PlaneWall>& getRightWall() const {return walls[1];}
    const std::shared_ptr<PlaneWall>& getFrontWall() const {return walls[3];}
    const std::shared_ptr<PlaneWall>& getBackWall() const {return walls[4];}

private:
    std::vector<std::shared_ptr<PlaneWall>> walls;
    Eigen::Vector3d lowerCorner;
    Eigen::Vector3d dimensions;

    // Helper methods to create walls
    void createWalls(const Eigen::Vector3d& velocity);
};


