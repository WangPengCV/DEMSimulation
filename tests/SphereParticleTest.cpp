#include "SphereParticle.h"
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
        auto manager = std::make_shared<ParticlePropertyManager>();
        PropertyTypeID id(1, 1);

        // Adding properties to manager
        manager->addSphereProperties(id, properties);

        SphereParticle SP(0,id,1,manager,Eigen::Vector3d(0,0,0));

        // Printing some initial properties of the sphere particle
        std::cout << "Sphere Particle Position: " << SP.getPosition().transpose() << std::endl;

        
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
