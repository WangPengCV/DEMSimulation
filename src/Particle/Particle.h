#pragma once
#include <Eigen/Dense>
#include "ParticlePropertyManager.h"
#include "PlaneWall.h"
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

    virtual double computeOverlap(const std::shared_ptr<Particle>& particle) = 0;

    virtual double computeOverlap(const std::shared_ptr<PlaneWall>& planewall) = 0;

    virtual void updateVelocity(double deltaTime) = 0;

    virtual void updateOmega(double deltaTime) = 0;

    virtual std::string save_tostring() const = 0;  

    void setPosition(Eigen::Vector3d &position);

    void setId(int id);


    void addForce(const Eigen::Vector3d &additionalForce);

    void resetForce();

    void addTorque(const Eigen::Vector3d &additionalTorque);

    void resetTorque();

    void updatePosition(double deltaTime);

    // Accessor methods (if needed)
    int getId() const { return id; }    
    const Eigen::Vector3d& getPosition() const { return position; }
    const PropertyTypeID& getType() const {return type;}
    std::shared_ptr<ParticlePropertyManager>  getParticlePropertyManager() const { return manager;}


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
