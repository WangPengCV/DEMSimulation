#include "ParticlePropertyManager.h"



void ParticlePropertyManager::addSphereProperties(const PropertyTypeID& id, std::shared_ptr<SphereProperties> properties)
{
    propertiesMap[id] = properties;
}
