#pragma once
#include "ParticleProperties.h"
<<<<<<< HEAD
=======
#include <cmath>
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200

class SphereProperties : public ParticleProperties {
public:
    // Use the constructor of the base class
    SphereProperties(double density, double mass, double radius, double rollingFriction, double slidingFriction,
                     double youngModulus, double restitution, double poissonRatio,double moment_of_inertia);
<<<<<<< HEAD


    void setMass(double Mass);
    void setRadius(double Radiys);
    void setMomentofinertia(double Momentofinertia);

    // Getter for the moment of inertia
    double getMomentOfInertia() const {
        return momentofinertia;
    }
    double getMass() const { return mass; }
    double getRadius() const { return radius; }
=======
    // Getter for the moment of inertia
    double getMomentOfInertia() const {
        return moment_of_inertia;
    }
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    std::string save_tostring() const override ;

private:

<<<<<<< HEAD
    double momentofinertia;
    double radius;
    double mass;

=======
    double moment_of_inertia;
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200

};
