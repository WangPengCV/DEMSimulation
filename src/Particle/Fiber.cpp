#include "Fiber.h"

Fiber::Fiber(int id, PropertyTypeID type, int startstate, int endstate, std::shared_ptr<ParticlePropertyManager>& manager)
: id(id),type(type),startstate(startstate),endstate(endstate),manager(manager)
{

}

void Fiber::setSphereCylinderBond(const std::vector<int>& spherecylinderbond)
{
    this->spherecylinderbond = spherecylinderbond;
}
const std::vector<int>& Fiber::getSphereCylinderBond() const
{
    return spherecylinderbond;
}
