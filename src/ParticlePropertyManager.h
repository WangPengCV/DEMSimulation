#pragma once
#include "ParticleProperties.h"
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
    void addProperties(const PropertyTypeID& id, std::shared_ptr<ParticleProperties> properties); 
   

    std::shared_ptr<ParticleProperties> getProperties(const PropertyTypeID& id) const 
    {
        auto it = propertiesMap.find(id);
        if (it != propertiesMap.end()) {
            return it->second;
        }
        return nullptr; // or throw an exception if preferred
    }

private:
    std::map<PropertyTypeID, std::shared_ptr<ParticleProperties>> propertiesMap;
};
