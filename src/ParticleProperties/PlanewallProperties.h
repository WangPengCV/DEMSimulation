#pragma once
#include "ParticleProperties.h"
#include <cmath>

class PlanewallProperties : public ParticleProperties {
public:
    // Use the constructor of the base class
<<<<<<< HEAD
    PlanewallProperties(double density, double thickness,  double rollingFriction, double slidingFriction,
=======
    PlanewallProperties(double density, double thickness, double mass, double rollingFriction, double slidingFriction,
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
                     double youngModulus, double restitution, double poissonRatio);
    virtual std::string save_tostring() const override;
    void setThickness(double thickness );
    double getThickness() const { return thickness; }
    
private:
    double thickness;
};
