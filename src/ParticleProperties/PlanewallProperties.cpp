#include "PlanewallProperties.h"

<<<<<<< HEAD
PlanewallProperties::PlanewallProperties(double density, double thickness, double rollingFriction, double slidingFriction,
                                         double youngModulus, double restitution, double poissonRatio)
    : ParticleProperties(density, rollingFriction, slidingFriction,
=======
PlanewallProperties::PlanewallProperties(double density, double thickness, double mass, double rollingFriction, double slidingFriction,
                                         double youngModulus, double restitution, double poissonRatio)
    : ParticleProperties(density, mass, thickness, rollingFriction, slidingFriction,
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
                         youngModulus, restitution, poissonRatio),
      thickness(thickness)
{
}

void PlanewallProperties::setThickness(double thickness)
{
<<<<<<< HEAD
    this->thickness = thickness;
=======
    thickness = thickness;
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
}

std::string PlanewallProperties::save_tostring() const
{
    std::ostringstream ss;
    ss << density << ", " << thickness << ", "
       << rolling_friction_coefficient << ", " << slide_friction_coefficient << ", " << Young_modulus << ", "
       << restitution << ", " << poisson_ratio;
    return ss.str();
}
