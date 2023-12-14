#pragma once
#include <set>
#include <tuple>
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