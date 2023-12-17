#include "SphereParticle.h"

SphereParticle::SphereParticle(int id, PropertyTypeID type, int state,
                               std::shared_ptr<ParticlePropertyManager> manager,
                               const Eigen::Vector3d &position,
                               const Eigen::Vector3d &velocity, const Eigen::Vector3d &omega,
                               const Eigen::Vector3d &force, const Eigen::Vector3d &torque)
    : Particle(id, type, state, manager, position, velocity, omega, force, torque)
{
}

SphereParticle::~SphereParticle()
{
    // Custom cleanup for SphereParticle, if needed
}

void SphereParticle::updateVelocity(double deltaTime)
{
    // Implement sphere-specific velocity update logic
    double mass = manager->getSphereProperties(type)->getMass();
    Eigen::Vector3d acceleration = force / mass;
    velocity += acceleration * deltaTime;
}

void SphereParticle::updateOmega(double deltaTime)
{
    // Implement sphere-specific angular velocity update logic
    double moment_of_inertia = manager->getSphereProperties(type)->getMomentOfInertia();
    Eigen::Vector3d angular_acceleration = torque / moment_of_inertia;
    omega += angular_acceleration * deltaTime;
}
