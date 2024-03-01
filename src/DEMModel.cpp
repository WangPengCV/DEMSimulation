#include "DEMModel.h"
#include <filesystem>

// Initialize DEMModel
DEMModel::DEMModel(const std::string &filename)
{
    DEMproperties = std::make_shared<DEMProperties>();
    DEMproperties->loadFromFile(filename);
    vis = std::make_shared<Visualization>(*DEMproperties);
    DEMproperties->initialSimulation();

    
}

void DEMModel::runSimulation()
{
    

    double currentTime = 0.0;
    double timeStep = DEMproperties->getTimestep();
    double totalTime = DEMproperties->getTotalTime();
    int showInterval = DEMproperties->getShowInterval();
    //generate particles
    int generate_number = 0;
    while (!DEMproperties->isGenerateComplete())
    {
        //DEMproperties->applyExternalForces();
        DEMproperties->handleCollisions();
        DEMproperties->bondforce();
        DEMproperties->motion();
        if(generate_number % 50000 == 0)
        {
            DEMproperties->generateRemainingParticles();
            vis->Update();
        }
        generate_number++;
    }
    //reach a quasi-static state
    generate_number =0;
    while(DEMproperties->getAverageVelocity() > 0.001 || generate_number < 100000)
    {
        //DEMproperties->applyExternalForces();
        DEMproperties->handleCollisions();
        DEMproperties->bondforce();
        DEMproperties->motion();
        if(generate_number % 50000 == 0)
        {
            vis->Update();
            std::cout << " average velocity is " << DEMproperties->getAverageVelocity() << std::endl;
        }
        generate_number++;
    }

    // task such as compression, damper......
    DEMproperties->initial_task();
    int iter_num = 0;
    std::string file_folder = "DEMProperties";
    std::filesystem::path dirPath(file_folder); // Directory path

    // Check if the directory exists
    if (!std::filesystem::exists(dirPath)) {
        // Create the directory if it does not exist
        std::filesystem::create_directories(dirPath);
    }
    while (currentTime < totalTime)
    {
        DEMproperties->applyExternalForces();
        DEMproperties->handleCollisions();
        DEMproperties->bondforce();

        DEMproperties->dotask();

        // Update visualization at specified intervals
        if (iter_num % showInterval == 0) {
            vis->Update();
            std::string filename = file_folder + "/DEMProperties" + std::to_string(iter_num) + ".dat";
            DEMproperties->saveToFile(filename);

        }
        currentTime += timeStep;
        iter_num++;
        DEMproperties->motion();



    }
}




