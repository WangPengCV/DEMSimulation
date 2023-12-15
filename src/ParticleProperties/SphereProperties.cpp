#pragma once
#include "SphereProperties.h"

SphereProperties::SphereProperties(double density, double mass, double radius, double rollingFriction, double slidingFriction,
                                   double youngModulus, double restitution, double poissonRatio)
    : ParticleProperties(density, mass, radius, rollingFriction, slidingFriction,
                         youngModulus, restitution, poissonRatio)
{
}

// Additional sphere-specific properties or methods
double SphereProperties::getVolume() const
{
    return (4.0 / 3.0) * PI * std::pow(getRadius(), 3);
}

double SphereProperties::getSurfaceArea() const
{
    return 4.0 * PI * std::pow(getRadius(), 2);
}
