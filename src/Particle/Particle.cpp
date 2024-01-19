#include "Particle.h"

Particle::Particle(int id, PropertyTypeID type, int state,
<<<<<<< HEAD
                   const std::shared_ptr<ParticlePropertyManager>& manager,
=======
                   std::shared_ptr<ParticlePropertyManager> manager,
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
                   const Eigen::Vector3d &position,
                   const Eigen::Vector3d &velocity, const Eigen::Vector3d &omega,
                   const Eigen::Vector3d &force, const Eigen::Vector3d &torque)
    : id(id), type(type), state(state), position(position), manager(manager),
      velocity(velocity), omega(omega), force(force), torque(torque)
{

}

Particle::~Particle()
{
    // Destructor implementation (if needed)
}

void Particle::setPosition(const Eigen::Vector3d &newPosition)
{
    position = newPosition;
}

void Particle::setId(int newId)
{
    id = newId;
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
