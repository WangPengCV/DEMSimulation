#pragma once
#include "SphereProperties.h"
#include <set>
#include <tuple>
#include <memory>
#include <map>

class PropertyTypeID {
public:
    PropertyTypeID(int category, int subType) : category(category), subType(subType) {}

    bool operator<(const PropertyTypeID& other) const {
        return std::tie(category, subType) < std::tie(other.category, other.subType);
    }

private:
    int category;  // For example, 1 for spheres, 2 for cylinders, etc.
    int subType;   // Subtype to distinguish different properties within the same category
};

class ParticlePropertyManager {
public:
   
    // Method to add SphereProperties
    void addSphereProperties(const PropertyTypeID& id, std::shared_ptr<SphereProperties> properties); 

    // Method to get SphereProperties
    std::shared_ptr<SphereProperties> getSphereProperties(const PropertyTypeID& id) const {
        auto it = propertiesMap.find(id);
        if (it != propertiesMap.end()) {
            // Use dynamic_cast to ensure the correct type is returned
            return std::dynamic_pointer_cast<SphereProperties>(it->second);
        }
        return nullptr;
    }
    

private:
    std::map<PropertyTypeID, std::shared_ptr<ParticleProperties>> propertiesMap;
};
