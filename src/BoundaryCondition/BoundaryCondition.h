#pragma once
#include <Eigen/Dense>
#include "ParticlePropertyManager.h"

class BoundaryCondition {
public:
    BoundaryCondition(int id, const PropertyTypeID& type, int state)
        : id(id), type(type), state(state) {}

    virtual ~BoundaryCondition() = default;

    // Accessors
    int getId() const { return id; }
    PropertyTypeID getType() const { return type; }
    int getState() const { return state; }

protected:
    int id;
    PropertyTypeID type;
    int state;
};
