#include "ParticlePropertyManager.h"

<<<<<<< HEAD
void ParticlePropertyManager::addSphereProperties(const PropertyTypeID& id, std::shared_ptr<SphereProperties>& properties)
{
    // Check if the property already exists
    if (propertiesMap.find(id) != propertiesMap.end()) {
        throw std::runtime_error("Property with the same ID already exists.");
    }    
    // Add the property if no conflict
    propertiesMap[id] = properties;
}

void ParticlePropertyManager::addPlanewallProperties(const PropertyTypeID& id, std::shared_ptr<PlanewallProperties>& properties)
{
     // Check if the property already exists
    if (propertiesMap.find(id) != propertiesMap.end()) {
        throw std::runtime_error("Property with the same ID already exists.");
    } 
    propertiesMap[id] = properties;
}

void ParticlePropertyManager::addFiberProperties(const PropertyTypeID &id, std::shared_ptr<FiberProperties>& properties)
{

     // Check if the property already exists
    if (propertiesMap.find(id) != propertiesMap.end()) {
        throw std::runtime_error("Property with the same ID already exists.");
    } 
=======


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
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    propertiesMap[id] = properties;
}
