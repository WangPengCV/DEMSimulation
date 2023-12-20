#include "DEMProperties.h"
DEMProperties::DEMProperties()
{
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
            line_process(line);
            parseLine(line);
        }
    }
}
bool DEMProperties::parseVector3d(std::istringstream &iss, Eigen::Vector3d &vec)
{
    for (int i = 0; i < 3; ++i)
    {
        if (!(iss >> vec[i]))
        {
            return false;
        }
    }
    return true;
}
void DEMProperties::parseLine(const std::string &line)
{
    // Parse the line and set properties accordingly
    // Example: "GRAVITY 0 -9.81 0" sets the gravity vector
    // Example: "PARTICLE SPHERE 1 0 0 0" defines a sphere particle
    // and so on...
    if (line.empty() || line[0] == '#')
        return; // Skip empty lines and comments

    std::istringstream iss(line);
    std::string entryType, token;
    iss >> entryType;

    if (entryType == "SPHERE_PROPERTIES")
    {
        parseSphereProperties(iss);
    }
    else if (entryType == "CYLINDER_PROPERTIES")
    {
    }
    else if (entryType == "PLANEWALL_PROPERTIES")
    {
        parsePlanewallProperties(iss);
    }
    else if (entryType == "BOUNDARY")
    {
        std::string boundaryType;
        iss >> boundaryType;

        if (boundaryType == "PLANEWALL")
        {
            int id, categoryTypeID, subTypeID, state;
            Eigen::Vector3d normal, corner1, corner2, corner3;
            parsePlaneWall(iss);
        }
        // Handle other boundary types similarly
    }
    else if (entryType == "RANDOM_PARTICLE")
    {
        parseRandomParticle(iss);
    }
    else if (entryType == "PARTICLE")
    {
        parseSpecificParticle(iss);
    }
    else if(entryType == "TIMESTEP")
    {
        iss >> timestep;
    }
    else if(entryType == "TOTAL_TIME")
    {
        iss >> totalTime;
    }
}
void DEMProperties::line_process(std::string &line)
{
    for (char &c : line)
    {
        if (c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n')
            c = ' ';
    }
    line.erase(0, line.find_first_not_of(" "));
    line.erase(line.find_last_not_of(" ") + 1);
}
void DEMProperties::parseSphereProperties(std::istringstream &iss)
{
    int category_id, subtype;
    double density, radius, rollingFriction, slidingFriction, youngModulus, restitution, poissonRatio;

    try
    {
        iss >> category_id >> subtype >> density >> radius >> rollingFriction >> slidingFriction >> youngModulus >> restitution >> poissonRatio;

        // Use the parsed values...
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing SPHERE_PROPERTIES: " << e.what() << std::endl;
        return;
    }
    // Calculations based on the properties
    double mass = 4.0 / 3.0 * PI * std::pow(radius, 3) * density;
    double moment_of_inertia = 2.0 / 5.0 * mass * std::pow(radius, 2);

    // Create SphereProperties object and store or process as needed
    auto sphereProperties = std::make_shared<SphereProperties>(density, mass, radius, rollingFriction, slidingFriction,
                                                               youngModulus, restitution, poissonRatio, moment_of_inertia);
    // Assuming you have a mechanism to add these properties to your manager
    manger->addSphereProperties(PropertyTypeID(category_id, subtype), sphereProperties);

    manger->addSphereType(category_id);
}
void DEMProperties::parsePlanewallProperties(std::istringstream &iss)
{
    int category_id, subtype;
    double density, thickness, rollingFriction, slidingFriction, youngModulus, restitution, poissonRatio;

    try
    {
        iss >> category_id >> subtype >> density >> thickness >> rollingFriction >> slidingFriction >> youngModulus >> restitution >> poissonRatio;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing Planewall_Propertie: " << e.what() << std::endl;
        return;
    }

    double mass = 0;

    auto properties = std::make_shared<PlanewallProperties>(density, thickness, mass, rollingFriction, slidingFriction,
                                                            youngModulus, restitution, poissonRatio);
    manger->addPlanewallType(category_id);
    manger->addPlanewallProperties(PropertyTypeID(category_id, subtype), properties);
}
void DEMProperties::parsePlaneWall(std::istringstream &iss)
{
    int id, category_id, subtype, state;
    Eigen::Vector3d normal, corner1, corner2, corner3, velocity;
    double meshResolution;
    try
    {
        iss >> id >> category_id >> subtype >> state;
        if (!parseVector3d(iss, normal) || !parseVector3d(iss, corner1) || !parseVector3d(iss, corner2) || !parseVector3d(iss, corner3) || !parseVector3d(iss, velocity))
        {
            throw std::runtime_error("Error parsing PLANEWALL boundary vectors.");
        }
        meshResolution = 0;
        auto pw = std::make_shared<PlaneWall>(id, PropertyTypeID(category_id, subtype), state, normal, corner1, corner2, corner3, velocity, meshResolution);
        planewalls.push_back(pw);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing PlaneWall: " << e.what() << std::endl;
        return;
    }
}
void DEMProperties::saveToFile(const std::string &filename) const
{

    std::ofstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Unable to open file for saving: " + filename);
    }

    for (const auto &pair : manger->getParticleProperties())
    {
        const auto &id = pair.first;
        const auto &prop = pair.second;
        std::map<int, ParticleType> typeMapping = manger->gettypeMapping();

        if (typeMapping[id.getCategory()] == ParticleType::SPHERE)
        {
            file << "SPHERE_PROPERTIES, " << id.getCategory() << ", " << id.getSubType() << ", " << prop->save_tostring() << "\n";
        }
        else if (typeMapping[id.getCategory()] == ParticleType::PLANEWALL)
        {
            file << "PLANEWALL_PROPERTIES, " << id.getCategory() << ", " << id.getSubType() << ", " << prop->save_tostring() << "\n";
        }
        // Add other property types as needed
    }

    // Write boundary conditions
    for (const auto &boundary : planewalls)
    {
        file << "BOUNDARY, " << boundary->save_tostring() << "\n";
    }

    for (const auto &particle : particles)
    {
        file << particle->save_tostring() << std::endl;
    }
    file.close();
}
void DEMProperties::parseRandomParticle(std::istringstream &iss)
{

    std::string type;
    int categoryId, subTypeId, state, count;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    iss >> type >> categoryId >> subTypeId >> state >> count >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;

    int id_index = particles.size();

    if (type == "SPHERE")
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> disX(xmin, xmax);
        std::uniform_real_distribution<> disY(ymin, ymax);
        std::uniform_real_distribution<> disZ(zmin, zmax);
        // Generate 'count' random particles within the specified ranges
        for (int i = 0; i < count; ++i)
        {
            double x, y, z;
            bool validPosition = false;
            do
            {
                validPosition = true;
                x = disX(gen);
                y = disY(gen);
                z = disZ(gen);
                auto tempsp = std::make_shared<SphereParticle>(id_index, PropertyTypeID(categoryId, subTypeId), state, manger, Eigen::Vector3d(x, y, z));
                for (const auto &plane : planewalls)
                {
                    double overlap_pw = tempsp->computeOverlap(plane);
                    if (overlap_pw > 0)
                    {
                        validPosition = false;
                        break;
                    }
                }

                if (validPosition)
                {
                    for (const auto &existparticle : particles)
                    {
                        double overlap_pp = tempsp->computeOverlap(existparticle);

                        if (overlap_pp > 0)
                        {
                            validPosition = false;
                            break;
                        }
                    }
                }

            } while (!validPosition);
            auto sp = std::make_shared<SphereParticle>(id_index, PropertyTypeID(categoryId, subTypeId), state, manger, Eigen::Vector3d(x, y, z));

            particles.push_back(sp);
            id_index++;
        }
    }
}
void DEMProperties::parseSpecificParticle(std::istringstream &iss)
{
    std::string type;
    int id, categoryId, subTypeId, state;
    double x, y, z, vx, vy, vz;
    iss >> type >> id >> categoryId >> subTypeId >> state >> x >> y >> z >> vx >> vy >> vz;
    if (type == "SPHERE")
    {
        // Create and add the specific particle to your simulation
        auto sp = std::make_shared<SphereParticle>(id, PropertyTypeID(categoryId, subTypeId), state, manger, Eigen::Vector3d(x, y, z), Eigen::Vector3d(vx, vy, vz));

        particles.push_back(sp);
    }
}