#pragma once
#include "BoundaryCondition.h"
#include "PlaneWall.h"
<<<<<<< HEAD
#include "ParticlePropertyManager.h"
=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
#include <unordered_map>
#include <Eigen/Dense>

class RectangularContainer : public BoundaryCondition {
public:
<<<<<<< HEAD
    RectangularContainer(){};
=======
    RectangularContainer();
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    // Constructor
    RectangularContainer(int id, const PropertyTypeID& type, int state, 
                         const Eigen::Vector3d& lowerCorner, 
                         const Eigen::Vector3d& dimensions, 
<<<<<<< HEAD
                         double rotation,
                         const Eigen::Vector3d &Center,
=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
                         const Eigen::Vector3d& velocity = Eigen::Vector3d::Zero());

    void rotateContainer(double angleDegrees, const Eigen::Vector3d& axis);

<<<<<<< HEAD
    void translateToCenter(const Eigen::Vector3d& simulationDomain);


    void setMass(std::shared_ptr<ParticlePropertyManager> manager);

    void setVelocity(const Eigen::Vector3d& Velocity);

    void setCenter(const Eigen::Vector3d& Center);

    void addForce(const Eigen::Vector3d& Force);

    void resetForce();

=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
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
<<<<<<< HEAD
    double getMass() const{return mass;}
    const Eigen::Vector3d& getVelocity() const {return velocity;}
    const Eigen::Vector3d& getForce() const {return force;}
    const Eigen::Vector3d& getCenter() const {return center;}

=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200

private:
    std::vector<std::shared_ptr<PlaneWall>> walls;
    Eigen::Vector3d lowerCorner;
    Eigen::Vector3d dimensions;
<<<<<<< HEAD
    Eigen::Vector3d center;

    double rotation;
    double mass;
    Eigen::Vector3d force;
    Eigen::Vector3d velocity;

    // Helper methods to create walls
    void createWalls();


=======

    // Helper methods to create walls
    void createWalls(const Eigen::Vector3d& velocity);
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
};


