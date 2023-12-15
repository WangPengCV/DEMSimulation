#include "ParticlePropertyManager.h"

void ParticlePropertyManager::addProperties(const PropertyTypeID &id, std::shared_ptr<ParticleProperties> properties)
{
    propertiesMap[id] = properties;

}

void ParticlePropertyManager::addSphereProperties(const PropertyTypeID& id, std::shared_ptr<SphereProperties> properties)
{
    propertiesMap[id] = properties;
}
