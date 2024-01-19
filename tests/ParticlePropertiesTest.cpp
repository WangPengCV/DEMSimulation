#include "ParticlePropertyManager.h" // Include the appropriate header
#include <iostream>

int main(int argc, char **argv)
{
    try
    {
        // Creating an instance of ParticleProperties
        auto properties = std::make_shared<SphereProperties>(7000, 0.5, 0.2, 0.01, 0.5, 10000000, 0.95, 0.3,10);
        // Printing some properties
        std::cout << "Density: " << properties->getDensity() << std::endl;
        std::cout << "Mass: " << properties->getMass() << std::endl;

        // Creating an instance of ParticlePropertyManager
        ParticlePropertyManager manager;
        PropertyTypeID id(1, 1);

        // Adding properties to manager
        manager.addSphereProperties(id, properties);

        // Retrieving and printing properties from manager
        auto retrievedProperties = manager.getSphereProperties(id);
        if (retrievedProperties)
        {
            std::cout << "Retrieved Mass: " << retrievedProperties->getMass() << std::endl;
        }
        else
        {
            std::cout << "Properties not found for given ID." << std::endl;
        }

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
