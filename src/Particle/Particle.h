#pragma once
#include <Eigen/Dense>
#include "ParticlePropertyManager.h"
#include <memory>

class Particle
{
public:
    Particle(int id, PropertyTypeID type, int state,
             std::shared_ptr<ParticlePropertyManager> manager,
             const Eigen::Vector3d &position,
             const Eigen::Vector3d &velocity = Eigen::Vector3d::Zero(),
             const Eigen::Vector3d &omega = Eigen::Vector3d::Zero(),
             const Eigen::Vector3d &force = Eigen::Vector3d::Zero(),
             const Eigen::Vector3d &torque = Eigen::Vector3d::Zero());

    virtual ~Particle();

    void addForce(const Eigen::Vector3d &additionalForce);
    void resetForce();
    void addTorque(const Eigen::Vector3d &additionalTorque);
    void resetTorque();
    void updatePosition(double deltaTime);
    virtual void updateVelocity(double deltaTime) = 0;
    virtual void updateOmega(double deltaTime) = 0;

    // Accessor methods (if needed)
    int getId() const { return id; }    
    const Eigen::Vector3d& getPosition() const { return position; }


    // ... other accessors ...

protected:
    int id;
    PropertyTypeID type;
    int state;

    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    Eigen::Vector3d omega;
    Eigen::Vector3d force;
    Eigen::Vector3d torque;

    std::shared_ptr<ParticlePropertyManager> manager;
};
