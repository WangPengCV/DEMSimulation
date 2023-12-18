#include "PlanewallProperties.h"

PlanewallProperties::PlanewallProperties(double density, double thickness, double mass, double rollingFriction, double slidingFriction,
                                         double youngModulus, double restitution, double poissonRatio)
    : ParticleProperties(density, mass, thickness, rollingFriction, slidingFriction,
                         youngModulus, restitution, poissonRatio),
      thickness(thickness)
{
}

void PlanewallProperties::setThickness(double thickness)
{
    thickness = thickness;
}

std::string PlanewallProperties::save_tostring() const
{
    std::ostringstream ss;
    ss << density << ", " << thickness << ", "
       << rolling_friction_coefficient << ", " << slide_friction_coefficient << ", " << Young_modulus << ", "
       << restitution << ", " << poisson_ratio;
    return ss.str();
}
