#pragma once
#include <Eigen/Dense>
#include "ParticlePropertyManager.h"
#include "PlaneWall.h"
#include <memory>

class Particle
{
public:
    Particle(int id, PropertyTypeID type, int state,
<<<<<<< HEAD
             const std::shared_ptr<ParticlePropertyManager>& manager,
=======
             std::shared_ptr<ParticlePropertyManager> manager,
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
             const Eigen::Vector3d &position,
             const Eigen::Vector3d &velocity = Eigen::Vector3d::Zero(),
             const Eigen::Vector3d &omega = Eigen::Vector3d::Zero(),
             const Eigen::Vector3d &force = Eigen::Vector3d::Zero(),
             const Eigen::Vector3d &torque = Eigen::Vector3d::Zero());
    
    virtual ~Particle();

<<<<<<< HEAD
=======
    virtual double computeOverlap(const std::shared_ptr<Particle>& particle) = 0;

    virtual double computeOverlap(const std::shared_ptr<PlaneWall>& planewall) = 0;

>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    virtual void updateVelocity(double deltaTime, Eigen::Vector3d& gravity) = 0;

    virtual void updateOmega(double deltaTime) = 0;

    virtual std::string save_tostring() const = 0;  

    void setPosition(const Eigen::Vector3d &position);

    void setId(int id);


    void addForce(const Eigen::Vector3d &additionalForce);

    void resetForce();

    void addTorque(const Eigen::Vector3d &additionalTorque);

    void resetTorque();

    void updatePosition(double deltaTime);

    // Accessor methods (if needed)
    int getId() const { return id; }    
    const Eigen::Vector3d& getPosition() const { return position; }
    const Eigen::Vector3d& getVelocity() const { return velocity; }
    const Eigen::Vector3d& getOmega() const { return omega; }


    const PropertyTypeID& getType() const {return type;}
    int getState() const {return state;}
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
