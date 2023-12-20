#include "DEMProperties.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include "Visualization.h"

int main(int argc, char **argv)
{
        std::string inputfilename = "InputFile.txt";
        DEMProperties demproperties;
       
        demproperties.loadFromFile(inputfilename);
        Visualization vis(demproperties);
        while(true)
        {
            vis.Update();  // Initial update to reflect the loaded data
        }

       
        std::shared_ptr<ParticlePropertyManager> manger = demproperties.getParticleManager();
        PropertyTypeID id(0, 0);

        // Retrieving and printing properties from manager
        auto retrievedProperties = manger->getSphereProperties(id);
        if (retrievedProperties)
        {
            std::cout << "Retrieved Mass: " << retrievedProperties->getMass() << std::endl;
        }
        else
        {
            std::cout << "Properties not found for given ID." << std::endl;
        }
        
        demproperties.saveToFile("OutputFile.txt");
        
    
    
        
    // Continue with rest of the program if file was loaded successfully
    return 0;
}