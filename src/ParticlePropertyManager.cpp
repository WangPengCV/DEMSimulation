#include "ParticlePropertyManager.h"

void ParticlePropertyManager::addProperties(const PropertyTypeID &id, std::shared_ptr<ParticleProperties> properties)
{
    propertiesMap[id] = properties;

}
