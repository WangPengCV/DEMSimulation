#pragma once
#include "ParticleProperties.h"
#include <cmath>

class SphereProperties : public ParticleProperties {
public:
    // Use the constructor of the base class
    SphereProperties(double density, double mass, double radius, double rollingFriction, double slidingFriction,
                     double youngModulus, double restitution, double poissonRatio);

    // Additional sphere-specific properties or methods
    double getVolume() const; 

    double getSurfaceArea() const;

};
