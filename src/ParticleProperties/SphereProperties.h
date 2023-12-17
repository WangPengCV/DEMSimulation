#pragma once
#include "ParticleProperties.h"
#include <cmath>

class SphereProperties : public ParticleProperties {
public:
    // Use the constructor of the base class
    SphereProperties(double density, double mass, double radius, double rollingFriction, double slidingFriction,
                     double youngModulus, double restitution, double poissonRatio,double moment_of_inertia);
                     // Getter for the moment of inertia
    double getMomentOfInertia() const {
        return moment_of_inertia;
    }
private:

    double moment_of_inertia;

};
