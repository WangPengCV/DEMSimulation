#pragma once
#include <Eigen/Dense>
class Particle 
{
public:
    Particle();
    virtual void step(double deltaTime) = 0; // Update particle state
protected:
};
