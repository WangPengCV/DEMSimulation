#include "DEMProperties.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

int main(int argc, char **argv)
{
    try
    {
        std::string inputfilename = "InputFile.txt";
        DEMProperties demproperties;

        demproperties.loadFromFile(inputfilename);

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
        
    }
    catch (const std::runtime_error &e)
    {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        // Handle the error, e.g., by exiting the program or asking the user for a different file
        return 1;
    }
    // Continue with rest of the program if file was loaded successfully
    return 0;
}