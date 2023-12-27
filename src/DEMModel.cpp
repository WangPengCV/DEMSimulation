#include "DEMModel.h"

// Initialize DEMModel
DEMModel::DEMModel(const std::string &filename)
{
    DEMproperties = std::make_shared<DEMProperties>();
    DEMproperties->loadFromFile(filename);
    vis = std::make_shared<Visualization>(*DEMproperties);


    
}

void DEMModel::runSimulation()
{
    double timeStep = DEMproperties->getTimestep();
    double totalTime = DEMproperties->getTotalTime();
    double currentTime = 0.0;
    DEMproperties->initialContactDetection();

    

    while (currentTime < totalTime)
    {
        DEMproperties->applyExternalForces();
        DEMproperties->handleCollisions();
        DEMproperties->motion();
        vis->Update();
        currentTime += timeStep;
    }
}




