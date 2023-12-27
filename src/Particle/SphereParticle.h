#pragma once
#include "Particle.h"

class SphereParticle : public Particle
{
public:
    SphereParticle(int id, PropertyTypeID type, int state,
                   std::shared_ptr<ParticlePropertyManager> manager,
                   const Eigen::Vector3d &position,
                   const Eigen::Vector3d &velocity = Eigen::Vector3d::Zero(),
                   const Eigen::Vector3d &omega = Eigen::Vector3d::Zero(),
                   const Eigen::Vector3d &force = Eigen::Vector3d::Zero(),
                   const Eigen::Vector3d &torque = Eigen::Vector3d::Zero());

    virtual ~SphereParticle();

    virtual double computeOverlap(const std::shared_ptr<Particle>& particle) override;
    virtual double computeOverlap(const std::shared_ptr<PlaneWall>& planewall) override;



    virtual std::string save_tostring() const override; 


    virtual void updateVelocity(double deltaTime, Eigen::Vector3d& gravity) override;
    virtual void updateOmega(double deltaTime) override;
};
