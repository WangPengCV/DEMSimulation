#pragma once
#include "SphereProperties.h"

SphereProperties::SphereProperties(double density, double mass, double radius, double rollingFriction, double slidingFriction,
                                   double youngModulus, double restitution, double poissonRatio,double moment_of_inertia)
    : ParticleProperties(density, mass, radius, rollingFriction, slidingFriction,
                         youngModulus, restitution, poissonRatio),moment_of_inertia(moment_of_inertia)
{
}


