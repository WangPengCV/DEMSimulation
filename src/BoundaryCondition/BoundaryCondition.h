#pragma once
#include <Eigen/Dense>
#include "ParticlePropertyManager.h"

class BoundaryCondition {
public:
    
    BoundaryCondition() = default;
    BoundaryCondition(int id, const PropertyTypeID& type, int state)
        : id(id), type(type), state(state) {}

    virtual ~BoundaryCondition() = default;
    virtual std::string save_tostring() const = 0;  // Pure virtual function

    void setState(int state);
    // Function to set movement type
    // Accessors
    int getId() const { return id; }
    PropertyTypeID getType() const { return type; }
    int getState() const { return state; }

protected:
    int id;
    PropertyTypeID type;
    int state;

};
