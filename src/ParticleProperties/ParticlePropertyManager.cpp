#include "ParticlePropertyManager.h"



void ParticlePropertyManager::addSphereProperties(const PropertyTypeID& id, std::shared_ptr<SphereProperties> properties)
{
    propertiesMap[id] = properties;
}
void ParticlePropertyManager::addSphereType(int category)
{
    typeMapping[category] = ParticleType::SPHERE;
}

void ParticlePropertyManager::addPlanewallType(int category)
{
    typeMapping[category] = ParticleType::PLANEWALL;
}


void ParticlePropertyManager::addPlanewallProperties(const PropertyTypeID& id, std::shared_ptr<PlanewallProperties> properties)
{
    propertiesMap[id] = properties;
}
