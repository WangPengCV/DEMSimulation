#include "DEMProperties.h"
DEMProperties::DEMProperties() {
    manger = std::make_shared<ParticlePropertyManager>();
}
void DEMProperties::loadFromFile(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Unable to open file: " + filename);
    }
    else
    {
        std::string line;
        while (std::getline(file, line))
        {
            parseLine(line);
        }
    }
}

void DEMProperties::parseLine(const std::string &line)
{
    // Parse the line and set properties accordingly
    // Example: "GRAVITY 0 -9.81 0" sets the gravity vector
    // Example: "PARTICLE SPHERE 1 0 0 0" defines a sphere particle
    // and so on...
    // Skip comment lines
    if (line[0] != '#')
    {
        std::istringstream iss(line);
        std::string entryType, token;
        getline(iss, entryType, ',');

        if (entryType == "SPHERE_PROPERTIES")
        {
            int category_id, subtype;
            double density, radius, rollingFriction, slidingFriction, youngModulus, restitution, poissonRatio;

            // Extract the data using getline with ',' as the delimiter
            getline(iss, token, ',');
            category_id = std::stoi(token);
            getline(iss, token, ',');
            subtype = std::stoi(token);
            getline(iss, token, ',');
            density = std::stod(token);
            getline(iss, token, ',');
            radius = std::stod(token);
            getline(iss, token, ',');
            rollingFriction = std::stod(token);
            getline(iss, token, ',');
            slidingFriction = std::stod(token);
            getline(iss, token, ',');
            youngModulus = std::stod(token);
            getline(iss, token, ',');
            restitution = std::stod(token);
            getline(iss, token);
            poissonRatio = std::stod(token); // No delimiter at the end

            // Check if all the properties were read successfully
            if (iss.fail())
            {
                std::cerr << "Failed to parse SPHERE_PROPERTIES." << std::endl;
                return; // Skip this line if there was an error
            }

            PropertyTypeID id(category_id, subtype);

            double mass = 4.0 / 3.0 * PI * radius * radius * radius * density;

            double moment_of_inertia = 2.0 / 5.0 * PI * radius * radius;

            auto properties = std::make_shared<SphereProperties>(density, mass, radius, rollingFriction, slidingFriction,
                                                                   youngModulus, restitution, poissonRatio,moment_of_inertia);

            manger->addSphereProperties(id, properties);
        }
        else if (entryType == "CYLINDER_PROPERTIES")
        {

            // ... Similarly, parse other entries
        }

    }
}
