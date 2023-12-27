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

void SphereParticle::updateVelocity(double deltaTime, Eigen::Vector3d& gravity)
{
    // Implement sphere-specific velocity update logic
    double mass = manager->getSphereProperties(type)->getMass();
    Eigen::Vector3d acceleration = force / mass + gravity;
    velocity += acceleration * deltaTime;
}

void SphereParticle::updateOmega(double deltaTime)
{
    // Implement sphere-specific angular velocity update logic
    double moment_of_inertia = manager->getSphereProperties(type)->getMomentOfInertia();
    Eigen::Vector3d angular_acceleration = torque / moment_of_inertia;
    omega += angular_acceleration * deltaTime;
}

double SphereParticle::computeOverlap(const std::shared_ptr<Particle> &another)
{
    int another_category = another->getType().getCategory();
    auto another_manager = another->getParticlePropertyManager();
    auto typemapping = another_manager->gettypeMapping();

    double overlap = 0;
    if (typemapping[another_category] == ParticleType::SPHERE)
    {
        auto radius1 = manager->getSphereProperties(this->getType())->getRadius();
        auto radius2 = another_manager->getSphereProperties(another->getType())->getRadius();

        // Compute the distance between the centers of the two spheres
        Eigen::Vector3d position1 = this->getPosition();
        Eigen::Vector3d position2 = another->getPosition();
        double distance = (position1 - position2).norm();

        // Compute the overlap (positive if spheres intersect, zero or negative otherwise)
        overlap = (radius1 + radius2) - distance;
       
    }
    return overlap > 0 ? overlap : 0.0;
}

double SphereParticle::computeOverlap(const std::shared_ptr<PlaneWall> &planewall)
{
    // Retrieve the sphere's center and radius
    Eigen::Vector3d sphereCenter = this->getPosition();
    double sphereRadius = manager->getSphereProperties(this->getType())->getRadius();

    // Retrieve the plane wall's point and normal
    Eigen::Vector3d planePoint = planewall->getCorner1();
    
    Eigen::Vector3d planeNormal = planewall->getNormal();

    // Compute the vector from a point on the plane to the sphere's center
    Eigen::Vector3d vecToSphereCenter = sphereCenter - planePoint;

    // Compute the distance from the sphere's center to the plane
    double distanceToPlane = vecToSphereCenter.dot(planeNormal);

    // Calculate the overlap (penetration depth)
    double overlap = sphereRadius - std::abs(distanceToPlane);

    return (overlap > 0.0) ? overlap : 0.0;
}



std::string SphereParticle::save_tostring() const {
    std::ostringstream ss;
    ss << "PARTICLE, "  << "SPHERE, "<< id << ", " << type.getCategory() << ", " << type.getSubType() << ", " << state << ", "
       << position.x() << ", " << position.y() << ", " << position.z() << ", "
       << velocity.x() << ", " << velocity.y() << ", " << velocity.z() << ", "
       << omega.x() << ", " << omega.y() << ", " << omega.z() ;
    return ss.str();
}