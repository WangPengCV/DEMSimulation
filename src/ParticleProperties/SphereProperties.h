#pragma once
#include "ParticleProperties.h"

class SphereProperties : public ParticleProperties {
public:
    // Use the constructor of the base class
    SphereProperties(double density, double mass, double radius, double rollingFriction, double slidingFriction,
                     double youngModulus, double restitution, double poissonRatio,double moment_of_inertia);


    void setMass(double Mass);
    void setRadius(double Radiys);
    void setMomentofinertia(double Momentofinertia);

    // Getter for the moment of inertia
    double getMomentOfInertia() const {
        return momentofinertia;
    }
    double getMass() const { return mass; }
    double getRadius() const { return radius; }
    std::string save_tostring() const override ;

private:

    double momentofinertia;
    double radius;
    double mass;


};
