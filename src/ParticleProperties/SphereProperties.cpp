<<<<<<< HEAD
=======
#pragma once
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
#include "SphereProperties.h"

SphereProperties::SphereProperties(double density, double mass, double radius, double rollingFriction, double slidingFriction,
                                   double youngModulus, double restitution, double poissonRatio, double moment_of_inertia)
<<<<<<< HEAD
    : ParticleProperties(density,  rollingFriction, slidingFriction,
                         youngModulus, restitution, poissonRatio),
                         radius(radius),mass(mass), momentofinertia(moment_of_inertia)
{
}

void SphereProperties::setMass(double Mass)
{
    mass = Mass;
    validateProperty(mass, "Mass");
}

void SphereProperties::setRadius(double Radius)
{
    radius = Radius;
    validateProperty(radius, "Radius");
}

void SphereProperties::setMomentofinertia(double Momentofinertia)
{
    momentofinertia = Momentofinertia;
    validateProperty(radius, "Radius");
}

=======
    : ParticleProperties(density, mass, radius, rollingFriction, slidingFriction,
                         youngModulus, restitution, poissonRatio),
      moment_of_inertia(moment_of_inertia)
{
}
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
std::string SphereProperties::save_tostring() const
{
    std::ostringstream ss;
    ss << density << ", " << radius << ", " << rolling_friction_coefficient << ", "
       << slide_friction_coefficient << ", " << Young_modulus << ", "
       << restitution << ", " << poisson_ratio;
    return ss.str();
}
