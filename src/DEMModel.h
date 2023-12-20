#pragma once
#include "Particle.h"
#include "SphereParticle.h"
#include "DEMProperties.h"
#include <vector>
#include <memory>

class DEMModel {
public:
    DEMModel();

    void initializeSimulation(const std::string& propertiesFile);
    
    void runSimulation(double totalTime, double timeStep);
    
    // Additional methods for managing particles, interactions, etc.
    void pauseSimulation();
    void resumeSimulation();
    void stopSimulation();

    void logSimulationData();
    void outputSimulationResults(const std::string& filename);

    void applyExternalForces();
    void addParticle(const Particle& particle);
    void removeParticle(int particleId);
private:

    std::shared_ptr<DEMProperties> DEMproperties;

    std::vector<std::shared_ptr<Particle>> particles;
    
    void updateParticleStates(double deltaTime);
    void handleCollisions();
    void handleBoundaryConditions();

    // Other private methods and member variables as needed
};
