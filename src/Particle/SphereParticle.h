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

    virtual void updateVelocity(double deltaTime) override;
    virtual void updateOmega(double deltaTime) override;
};
