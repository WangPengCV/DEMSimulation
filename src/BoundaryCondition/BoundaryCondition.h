#pragma once
#include <Eigen/Dense>
#include "ParticlePropertyManager.h"

class BoundaryCondition {
public:
<<<<<<< HEAD
    
    BoundaryCondition() = default;
=======
    BoundaryCondition(){}
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    BoundaryCondition(int id, const PropertyTypeID& type, int state)
        : id(id), type(type), state(state) {}

    virtual ~BoundaryCondition() = default;
    virtual std::string save_tostring() const = 0;  // Pure virtual function

    void setState(int state);
<<<<<<< HEAD
    // Function to set movement type
=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    // Accessors
    int getId() const { return id; }
    PropertyTypeID getType() const { return type; }
    int getState() const { return state; }

protected:
    int id;
    PropertyTypeID type;
    int state;
<<<<<<< HEAD

=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
};
