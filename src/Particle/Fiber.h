#pragma once
#include "ParticlePropertyManager.h"

class Fiber
{
public:
    Fiber(int id, PropertyTypeID type, int startstate, int endstate, std::shared_ptr<ParticlePropertyManager>& manager);
    const std::vector<int>& getSphereCylinderBond() const;
    void setSphereCylinderBond( const std::vector<int>& spherecylinderbond);
private:
    int id;
    PropertyTypeID type;
    int startstate;
    int endstate;
    std::shared_ptr<ParticlePropertyManager> manager;
    std::vector<int> spherecylinderbond;

};