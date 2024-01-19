#include "ParticleProperties.h"

<<<<<<< HEAD
ParticleProperties::ParticleProperties(double density, double rollingFriction, double slidingFriction,
                                       double youngModulus, double restitution, double poissonRatio)
    : density(density), rolling_friction_coefficient(rollingFriction),
=======
ParticleProperties::ParticleProperties(double density, double mass, double radius, double rollingFriction, double slidingFriction,
                                       double youngModulus, double restitution, double poissonRatio)
    : density(density), mass(mass), radius(radius), rolling_friction_coefficient(rollingFriction),
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
      slide_friction_coefficient(slidingFriction), Young_modulus(youngModulus),
      restitution(restitution), poisson_ratio(poissonRatio)
{
    validateProperties();
}

void ParticleProperties::setDensity(double Density)
{
    density = Density;
    validateProperty(density, "Density");
}

<<<<<<< HEAD

=======
void ParticleProperties::setMass(double Mass)
{
    mass = Mass;
    validateProperty(mass, "Mass");
}

void ParticleProperties::setRadius(double Radius)
{
    radius = Radius;
    validateProperty(radius, "Radius");
}
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200

void ParticleProperties::setRollingFriction(double RollingFriction)
{
    rolling_friction_coefficient = RollingFriction;
    validateProperty(rolling_friction_coefficient, "RollingFriction");
}

void ParticleProperties::setSlidingFriction(double SlidingFriction)
{
    slide_friction_coefficient = SlidingFriction;
    validateProperty(slide_friction_coefficient, "SlidingFriction");

}

void ParticleProperties::setYoungModulus(double YoungModulus)
{
    Young_modulus = YoungModulus;
    validateProperty(Young_modulus, "YoungModulus");

}

void ParticleProperties::setRestitution(double Restitution)
{
    restitution = Restitution;
    validateProperty(restitution, "Restitution");

}

void ParticleProperties::setPoissonRatio(double PoissonRatio)
{
    poisson_ratio = PoissonRatio;
    validateProperty(poisson_ratio, "PoissonRatio");

}

void ParticleProperties::validateProperties()
{
    validateProperty(density, "Density");
<<<<<<< HEAD
=======
    validateProperty(mass, "Mass");
    validateProperty(radius, "Radius");

>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    validateProperty(rolling_friction_coefficient, "RollingFriction");
    validateProperty(slide_friction_coefficient, "SlidingFriction");
    validateProperty(Young_modulus, "YoungModulus");
    validateProperty(restitution, "Restitution");
    validateProperty(poisson_ratio, "PoissonRatio");
}

void ParticleProperties::validateProperty(double value, const std::string &name)
{
<<<<<<< HEAD
     if (value < 0) {
        throw std::invalid_argument(name + " must be positive. Given: " + std::to_string(value));
=======
    if (value < 0)
    {
        throw std::invalid_argument(name + " must be positive.");
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    }
}
