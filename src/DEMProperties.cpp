#include "DEMProperties.h"
DEMProperties::DEMProperties()
{
    manager = std::make_shared<ParticlePropertyManager>();
    contactforce = std::make_shared<ContactForce>();
    gridbasedcontactdetection = std::make_shared<GridBasedContactDetection>();
    globalforces = Eigen::Vector3d::Zero();
    gravity = Eigen::Vector3d::Zero();
    gernerateFlag = true;
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
            parsePlaneWall(iss);
        }
        else if(boundaryType == "RECTANGULARCONTAINER")
        {
            parseRectangularContainer(iss);
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
    else if (entryType == "TIMESTEP")
    {
        iss >> timestep;
    }
    else if (entryType == "TOTAL_TIME")
    {
        iss >> totalTime;
    }
    else if (entryType == "GRAVITY")
    {
        iss >> gravity.x() >> gravity.y() >> gravity.z();
    }
    else if (entryType == "FORCE")
    {
        std::string forceType;
        iss >> forceType;
        if (forceType == "SPECIFIC")
        {
            Eigen::Vector3d force;
            int particle_id;
            iss >> particle_id >> force.x() >> force.y() >> force.z();
            specificforces[particle_id] = force;
        }
        else if (forceType == "GLOBAL")
        {
            iss >> globalforces.x() >> globalforces.y() >> globalforces.z();
        }
    }
    else if (entryType == "DIMENSIONS")
    {
        iss >> simulationdimensions.x() >> simulationdimensions.y() >> simulationdimensions.z();
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
    manager->addSphereProperties(PropertyTypeID(category_id, subtype), sphereProperties);

    manager->addSphereType(category_id);
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
    manager->addPlanewallType(category_id);
    manager->addPlanewallProperties(PropertyTypeID(category_id, subtype), properties);
}
void DEMProperties::parsePlaneWall(std::istringstream &iss)
{
    int id, category_id, subtype, state;
    Eigen::Vector3d normal, corner1, corner2, corner3, velocity;
    try
    {
        iss >> id >> category_id >> subtype >> state;
        if (!parseVector3d(iss, normal) || !parseVector3d(iss, corner1) || !parseVector3d(iss, corner2) || !parseVector3d(iss, corner3) || !parseVector3d(iss, velocity))
        {
            throw std::runtime_error("Error parsing PLANEWALL boundary vectors.");
        }
        auto pw = std::make_shared<PlaneWall>(id, PropertyTypeID(category_id, subtype), state, normal, corner1, corner2, corner3, velocity);
        planewalls.push_back(pw);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing PlaneWall: " << e.what() << std::endl;
        return;
    }
}
void DEMProperties::parseRectangularContainer(std::istringstream &iss)
{
    int id, category_id, subtype, state;
    Eigen::Vector3d lowerCorner, dimensions, velocity;
    try
    {
        iss >> id >> category_id >> subtype >> state;
        if (!parseVector3d(iss, lowerCorner) || !parseVector3d(iss, dimensions)  || !parseVector3d(iss, velocity))
        {
            throw std::runtime_error("Error parsing PLANEWALL boundary vectors.");
        }
        auto rc = std::make_shared<RectangularContainer>(id, PropertyTypeID(category_id, subtype), state,lowerCorner,dimensions,velocity);
        rc->rotateContainer(45,Eigen::Vector3d(0,0,1));
        rectangularcontainer = rc;

        planewalls.push_back(rectangularcontainer->getLeftWall());
        planewalls.push_back(rectangularcontainer->getRightWall());
        planewalls.push_back(rectangularcontainer->getBottomWall());
        planewalls.push_back(rectangularcontainer->getFrontWall());
        planewalls.push_back(rectangularcontainer->getBackWall());
        planewalls.push_back(rectangularcontainer->getTopWall());
       
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing RectangularContainer: " << e.what() << std::endl;
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
    file << "DIMENSIONS,  " << simulationdimensions.x() << ", " << simulationdimensions.y() << ", " << simulationdimensions.z() << "\n";
    file << "TIMESTEP, " << timestep << "\n";
    file << "TOTAL_TIME, " << totalTime << "\n";
    file << "GRAVITY, " << gravity.x() << ", " << gravity.y() << ", " << gravity.z() << "\n";

    for (const auto &pair : specificforces)
    {
        const auto &id = pair.first;
        const auto &force = pair.second;

        file << "FORCE, "
             << "SPECIFIC, " << id << "," << force.x() << ", " << force.y() << ", " << force.z() << "\n";
    }

    file << "FORCE, "
         << "GLOBAL, , " << globalforces.x() << ", " << globalforces.y() << ", " << globalforces.z() << "\n";

    for (const auto &pair : manager->getParticleProperties())
    {
        const auto &id = pair.first;
        const auto &prop = pair.second;
        std::map<int, ParticleType> typeMapping = manager->gettypeMapping();

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

   
    
    file << "BOUNDARY, " << rectangularcontainer->save_tostring() << "\n";
    

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
        gernerateFlag= true;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> disX(xmin, xmax);
        std::uniform_real_distribution<> disY(ymin, ymax);
        std::uniform_real_distribution<> disZ(zmin, zmax);
        int failedSpheres = 0;
        // Generate 'count' random particles within the specified ranges
        for (int i = 0; i < count; ++i)
        {
            double x, y, z;
            bool validPosition = false;
            int count_number = 0;
            do
            {
                count_number++;
                validPosition = true;
                x = disX(gen);
                y = disY(gen);
                z = disZ(gen);
                auto tempsp = std::make_shared<SphereParticle>(id_index, PropertyTypeID(categoryId, subTypeId), state, manager, Eigen::Vector3d(x, y, z));
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

            } while (!validPosition && count_number < 100);
            if (validPosition)
            {
                auto sp = std::make_shared<SphereParticle>(id_index, PropertyTypeID(categoryId, subTypeId), state, manager, Eigen::Vector3d(x, y, z));
                particles.push_back(sp);
                id_index++;
            }
            else
            {
                gernerateFlag = false;
                failedSpheres = count - i;
                // Use an ostringstream to format the data as a string
                std::ostringstream oss;
                oss << type << " "
                    << categoryId << " "
                    << subTypeId << " "
                    << state << " "
                    << failedSpheres << " "
                    << xmin << " "
                    << xmax << " "
                    << ymin << " "
                    << ymax << " "
                    << zmin << " "
                    << zmax;
                
                gernerateInfor[categoryId] = oss.str();

                break;
            }
        }
        // You can check generateFlag and failedParticles here for further actions
        if (!gernerateFlag)
        {
            std::cout << "Can't generate all spheres. Failed spheres: " << failedSpheres << std::endl;
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
        auto sp = std::make_shared<SphereParticle>(id, PropertyTypeID(categoryId, subTypeId), state, manager, Eigen::Vector3d(x, y, z), Eigen::Vector3d(vx, vy, vz));

        particles.push_back(sp);
    }
}
void DEMProperties::generateRemainingParticles()
{
    if(!gernerateFlag)
    {
        auto typemapping = manager->gettypeMapping();

        for(const auto& sub : gernerateInfor)
        {
            if(typemapping[sub.first] == ParticleType::SPHERE)
            {
                std::istringstream iss(sub.second);
                parseRandomParticle(iss);
            }
           
        }
    }
}
void DEMProperties::initialContactDetection()
{
    contactforce->addParticleProperties(manager);
    double minRadius = std::numeric_limits<double>::max();

    for (const auto &particleproperties : manager->getParticleProperties())
    {
        auto typeMapping = manager->gettypeMapping();
        if (typeMapping[particleproperties.first.getCategory()] == ParticleType::SPHERE)
        {
            double radius = particleproperties.second->getRadius();
            if (radius < minRadius)
            {
                minRadius = radius;
            }
        }
    }
    gridbasedcontactdetection->initial(simulationdimensions.x(), simulationdimensions.y(), simulationdimensions.z(), minRadius);
}
void DEMProperties::handleCollisions()
{
    std::unordered_map<int, std::unordered_set<int>> pp_contact_paris;
    std::unordered_map<int, std::unordered_set<int>> pw_contact_paris;



    gridbasedcontactdetection->ParticleBroadPhase(particles, pp_contact_paris);
    gridbasedcontactdetection->planewallBroadPhase(particles, planewalls, pw_contact_paris);

    auto typemapping = manager->gettypeMapping();

    for (const auto &contact_list : pp_contact_paris)
    {
        auto id1 = contact_list.first;
        PropertyTypeID type1 = particles[id1]->getType();
        for (auto &id2 : contact_list.second)
        {
            PropertyTypeID type2 = particles[id2]->getType();

            if (typemapping[type1.getCategory()] == ParticleType::SPHERE && typemapping[type2.getCategory()] == ParticleType::SPHERE)
            {
                auto sp1 = std::static_pointer_cast<SphereParticle>(particles[id1]);
                auto sp2 = std::static_pointer_cast<SphereParticle>(particles[id2]);

                contactforce->computeSphereSphereForce(sp1, sp2, timestep);
            }
        }
    }

    for (const auto &contact_list : pw_contact_paris)
    {
        auto wall_id = contact_list.first;
        for (auto &id : contact_list.second)
        {
            PropertyTypeID type = particles[id]->getType();

            if (typemapping[type.getCategory()] == ParticleType::SPHERE)
            {
                auto wall = planewalls[wall_id];
                auto sphere = std::static_pointer_cast<SphereParticle>(particles[id]);

                contactforce->computePlaneWallSphereForce(wall, sphere, timestep);
            }
        }
    }

   
}
void DEMProperties::applyExternalForces()
{
    for (auto &particle : particles)
    {

        particle->addForce(globalforces);
    }
    for (auto &indexforce : specificforces)
    {
        particles[indexforce.first]->addForce(indexforce.second);
    }
    // add gravity force
}
void DEMProperties::motion()
{
    for (auto &particle : particles)
    {
        if (particle->getState() == 1)
        {
            particle->updateVelocity(timestep, gravity);
            particle->updateOmega(timestep);
            particle->updatePosition(timestep);
        }
        particle->resetForce();
        particle->resetTorque();
    }
    for (auto &planewall : planewalls)
    {
        if (planewall->getState() == 1)
        {
            Eigen::Vector3d v = planewall->getVelociy();
            auto type = planewall->getType();
            v += planewall->getForce() / manager->getPlanewallProperties(type)->getMass() * timestep;
            Eigen::Vector3d p1 = planewall->getCorner1();
            Eigen::Vector3d p2 = planewall->getCorner2();
            Eigen::Vector3d p3 = planewall->getCorner3();
            Eigen::Vector3d displancement = v * timestep;
            p1 += displancement;
            p2 += displancement;
            p3 += displancement;
        }
        planewall->resetForce();
    }
    contactforce->updateContactInformation();
}