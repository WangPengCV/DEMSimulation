#pragma once
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <sstream>

const double PI = 3.1415926535897932384626;
class ParticleProperties
{
public:
<<<<<<< HEAD
    ParticleProperties(double density, double rollingFriction, double slidingFriction,
                       double youngModulus, double restitution, double poissonRatio);
    ParticleProperties(){}
=======
    ParticleProperties(double density, double mass, double radius, double rollingFriction, double slidingFriction,
                       double youngModulus, double restitution, double poissonRatio);
    ParticleProperties();
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    // Add a virtual destructor
    virtual ~ParticleProperties() = default;

    virtual std::string save_tostring() const = 0;  // Pure virtual function


    
    void setDensity(double density );
<<<<<<< HEAD
=======
    void setMass(double mass );
    void setRadius( double radius);

>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    void setRollingFriction(double rollingFriction);
    void setSlidingFriction(double slidingFriction);
    void setYoungModulus(double youngModulus);
    void setRestitution(double restitution);
    void setPoissonRatio(double poissonRatio);



    double getDensity() const { return density; }
<<<<<<< HEAD
=======
    double getMass() const { return mass; }
    double getRadius() const { return radius; }

>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    double getRollingFriction() const { return rolling_friction_coefficient; }
    double getSlidingFriction() const { return slide_friction_coefficient; }
    double getYoungModulus() const { return Young_modulus; }
    double getRestitution() const { return restitution; }
    double getPoissonRatio() const { return poisson_ratio; }

protected:
    void validateProperties();

    void validateProperty(double value, const std::string &name);

    // Geometric properties
    double density;
<<<<<<< HEAD
=======
    double mass;
    double radius;
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    double rolling_friction_coefficient;
    double slide_friction_coefficient;

    // Contact properties
    double Young_modulus;
    double restitution;
    double poisson_ratio;
};