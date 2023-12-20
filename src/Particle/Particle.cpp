#include "Particle.h"

Particle::Particle(int id, PropertyTypeID type, int state,
                   std::shared_ptr<ParticlePropertyManager> manager,
                   const Eigen::Vector3d &position,
                   const Eigen::Vector3d &velocity, const Eigen::Vector3d &omega,
                   const Eigen::Vector3d &force, const Eigen::Vector3d &torque)
    : id(id), type(type), state(state), position(position), manager(manager),
      velocity(velocity), omega(omega), force(force), torque(torque)
{
    manager = std::make_shared<ParticlePropertyManager>();
}

Particle::~Particle()
{
    // Destructor implementation (if needed)
}

void Particle::setPosition(Eigen::Vector3d &position)
{
    position = position;
}

void Particle::setId(int id)
{
    id = id;
}

void Particle::addForce(const Eigen::Vector3d &additionalForce)
{
    force += additionalForce;
}

void Particle::resetForce()
{
    force.setZero();
}

void Particle::addTorque(const Eigen::Vector3d &additionalTorque)
{
    torque += additionalTorque;
}

void Particle::resetTorque()
{
    torque.setZero();
}

void Particle::updatePosition(double deltaTime)
{
    position += velocity * deltaTime;
}
