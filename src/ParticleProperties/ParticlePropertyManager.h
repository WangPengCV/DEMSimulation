#pragma once
#include "SphereProperties.h"
#include "PlanewallProperties.h"
#include <set>
#include <tuple>
#include <memory>
#include <map>
enum class ParticleType {
    SPHERE,
    PLANEWALL
};
class PropertyTypeID {
public:
    PropertyTypeID(){}  
    PropertyTypeID(int category, int subType) : category(category), subType(subType) {}

    bool operator<(const PropertyTypeID& other) const {
        return std::tie(category, subType) < std::tie(other.category, other.subType);
    }
    int getCategory() const { return category; }
    int getSubType() const {return subType;}


private:
    int category;  // For example, 1 for spheres, 2 for cylinders, etc.
    int subType;   // Subtype to distinguish different properties within the same category
};

class ParticlePropertyManager {
public:
   
    // Method to add SphereProperties
    void addSphereProperties(const PropertyTypeID& id, std::shared_ptr<SphereProperties> properties); 
    void addSphereType(int category);
    // Method to add PlanewallProperties
    void addPlanewallProperties(const PropertyTypeID& id, std::shared_ptr<PlanewallProperties> properties); 
    void addPlanewallType(int category);
    // Method to get SphereProperties
    std::shared_ptr<SphereProperties> getSphereProperties(const PropertyTypeID& id) const {
        auto it = propertiesMap.find(id);
        if (it != propertiesMap.end()) {
            // Use dynamic_cast to ensure the correct type is returned
            return std::dynamic_pointer_cast<SphereProperties>(it->second);
        }
        return nullptr;
    }

    // Method to get SphereProperties
    std::shared_ptr<PlanewallProperties> getPlanewallProperties(const PropertyTypeID& id) const {
        auto it = propertiesMap.find(id);
        if (it != propertiesMap.end()) {
            // Use dynamic_cast to ensure the correct type is returned
            return std::dynamic_pointer_cast<PlanewallProperties>(it->second);
        }
        return nullptr;
    }
    
    const std::map<PropertyTypeID, std::shared_ptr<ParticleProperties>>& getParticleProperties() const{
        return propertiesMap;
    }
    const std::map<int, ParticleType>& gettypeMapping() const{
        return typeMapping;
    }

private:
    std::map<PropertyTypeID, std::shared_ptr<ParticleProperties>> propertiesMap;
    std::map<int, ParticleType> typeMapping;
};
