#pragma once
#include "Particle.h"
#include "SphereParticle.h"
#include "DEMProperties.h"
#include "Visualization.h"
#include <vector>
#include <memory>

class DEMModel {
public:
    DEMModel(const std::string& filename);

    void runSimulation();
    
   
private:

    std::shared_ptr<DEMProperties> DEMproperties;
    std::shared_ptr<Visualization> vis;

};
