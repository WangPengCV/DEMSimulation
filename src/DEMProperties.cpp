#include "DEMProperties.h"
DEMProperties::DEMProperties()
{
    manager = std::make_shared<ParticlePropertyManager>();
    contactforce = std::make_shared<ContactForce>();
    gridbasedcontactdetection = std::make_shared<GridBasedContactDetection>();
    globalforces = Eigen::Vector3d::Zero();
    gravity = Eigen::Vector3d::Zero();
    gernerateSphereFlag = true;
    gernerateFiberFlag = true;
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
    else if (entryType == "FIBER_PROPERTIES")
    {
        parseFiberProperties(iss);
    }
    else if (entryType == "BOUNDARY")
    {
        std::string boundaryType;
        iss >> boundaryType;

        if (boundaryType == "PLANEWALL")
        {
            parsePlaneWall(iss);
        }
        else if (boundaryType == "RECTANGULARCONTAINER")
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
    else if (entryType == "SHOWINTERVAL")
    {
        iss >> showInterval;
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
    else if (entryType == "TASK")
    {
        std::string tasktype;
        iss >> tasktype;

        if (tasktype == "BOUNDARY")
        {
            std::string boundarytype, movetype;
            int id;
            iss >> boundarytype;
            if (boundarytype == "RECTANGULARCONTAINER")
            {
                iss >> id >> movetype;
                auto &subboundarytask = dynamicBoundaryTask[boundarytype];
                if (movetype == "SPRINGOSCILLATOR")
                {
                    double k, c, p;
                    iss >> k >> c >> p;
                    std::shared_ptr<DynamicBoundaryMover> dynamicMover = std::make_shared<SpringOscillator>(k, c, p);
                    subboundarytask.insert({id, dynamicMover});
                }
            }
            // double springk, springc, compression;
            // iss >> springk >> springc >> compression;
            // taskInfor[task].push_back(springk);
            // taskInfor[task].push_back(springc);
            // taskInfor[task].push_back(compression);
        }
        else if (tasktype == "RANDOM_PACKING")
        {
            // taskInfor[task].push_back(0);
        }
        else if (tasktype == "EVENT")
        {
            std::string particletype, movetype;
            int id;
            iss >> particletype;
            if (particletype == "SPHERE")
            {
                iss >> id >> movetype;

                if (movetype == "FREEMOTIONTASK")
                {
                    Eigen::Vector3d position, velocity;
                    parseVector3d(iss, position);
                    parseVector3d(iss, velocity);
                    auto &subeventtask = eventTasks[particletype];

                    std::shared_ptr<EventTask> dynamicMover = std::make_shared<FreeMotionTask>(position, velocity);
                    subeventtask.insert({id, dynamicMover});
                }
            }
        }
    }
    else if (entryType == "SphereToSphereContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();

        contactforce->getSphereSphereContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    else if (entryType == "SphereToPlanewallContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();

        contactforce->getWallSphereContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    else if (entryType == "SphereToFiberContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();
        contactforce->getSphereFiberContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    else if (entryType == "FiberToFiberContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();

        // Assuming a method to store fiber to fiber contact information
        contactforce->getFiberFiberContactInformationList()[info.particleId1][info.particleId2] = info;
    }
    else if (entryType == "WallToFiberContact")
    {
        ContactInformation info;
        iss >> info.particleId1 >> info.particleId2;
        iss >> info.normalForce.x() >> info.normalForce.y() >> info.normalForce.z();
        iss >> info.tangentialForce.x() >> info.tangentialForce.y() >> info.tangentialForce.z();
        iss >> info.tangentialDisplacement.x() >> info.tangentialDisplacement.y() >> info.tangentialDisplacement.z();
        iss >> info.previousTangentialVelocity.x() >> info.previousTangentialVelocity.y() >> info.previousTangentialVelocity.z();

        // Store the contact information
        contactforce->getWallFiberContactInformationList()[info.particleId1][info.particleId2] = info;
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

    auto properties = std::make_shared<PlanewallProperties>(density, thickness, rollingFriction, slidingFriction,
                                                            youngModulus, restitution, poissonRatio);
    manager->addPlanewallProperties(PropertyTypeID(category_id, subtype), properties);
}
void DEMProperties::parseFiberProperties(std::istringstream &iss)
{
    int category_id, subtype, nodenumber;
    double density, radius, rollingFriction, slidingFriction, youngModulus, restitution,
        poissonRatio, normalmodulus, shearmodulus, twistmodulus, bendingmodulus, aspectratio, bonddampingcoefficient;

    try
    {
        iss >> category_id >> subtype >> density >> radius >> rollingFriction >> slidingFriction >> youngModulus >> restitution >> poissonRatio >> normalmodulus >> shearmodulus >> twistmodulus >> bendingmodulus >> nodenumber >> aspectratio >> bonddampingcoefficient;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing Planewall_Propertie: " << e.what() << std::endl;
        return;
    }
    double elementlength = (aspectratio - 1) * 2 * radius / (nodenumber - 1);
    auto properties = std::make_shared<FiberProperties>(density, rollingFriction, slidingFriction,
                                                        youngModulus, restitution, poissonRatio, radius, elementlength, normalmodulus,
                                                        shearmodulus, twistmodulus, bendingmodulus, nodenumber, aspectratio, bonddampingcoefficient);
    manager->addFiberProperties(PropertyTypeID(category_id, subtype), properties);
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
    int id, category_id, subtype, state, rotate;
    Eigen::Vector3d lowerCorner, dimensions, velocity, center;
    try
    {
        iss >> id >> category_id >> subtype >> state;
        if (!parseVector3d(iss, lowerCorner) || !parseVector3d(iss, dimensions))
        {
            throw std::runtime_error("Error parsing PLANEWALL boundary vectors.");
        }
        iss >> rotate;
        parseVector3d(iss, center);
        parseVector3d(iss, velocity);
        auto rc = std::make_shared<RectangularContainer>(id, PropertyTypeID(category_id, subtype), state, lowerCorner, dimensions, rotate, center, velocity);
        rc->rotateContainer(rotate, Eigen::Vector3d(0, 0, 1));
        rc->translateToCenter(simulationdimensions);
        rectangularcontainer = rc;
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
    file << "SHOWINTERVAL, " << showInterval << "\n";
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

        if (auto sphereProerties = std::dynamic_pointer_cast<SphereProperties>(prop))
        {
            file << "SPHERE_PROPERTIES, " << id.getCategory() << ", " << id.getSubType() << ", " << prop->save_tostring() << "\n";
        }
        else if (auto planewallProerties = std::dynamic_pointer_cast<PlanewallProperties>(prop))
        {
            file << "PLANEWALL_PROPERTIES, " << id.getCategory() << ", " << id.getSubType() << ", " << prop->save_tostring() << "\n";
        }
        else if (auto planewallProerties = std::dynamic_pointer_cast<FiberProperties>(prop))
        {
            file << "FIBER_PROPERTIES, " << id.getCategory() << ", " << id.getSubType() << ", " << prop->save_tostring() << "\n";
        }
        // Add other property types as needed
    }

    // Write boundary conditions
    for (const auto &boundary : planewalls)
    {
        file << "BOUNDARY, " << boundary->save_tostring() << "\n";
    }

    if (rectangularcontainer)
    {
        file << "BOUNDARY, " << rectangularcontainer->save_tostring() << "\n";
    }

    for (const auto &sphere : sphereparticles)
    {
        file << sphere->save_tostring() << std::endl;
    }
    // Saving sphere-to-sphere contact information
    for (const auto &pair1 : contactforce->getSphereSphereContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "SphereToSphereContact"
                 << ", "
                 << info.particleId1 << ", "
                 << info.particleId2 << ", "
                 << info.normalForce.transpose() << ", "
                 << info.tangentialForce.transpose() << ", "
                 << info.tangentialDisplacement.transpose() << ", "
                 << info.previousTangentialVelocity.transpose() << "\n";
        }
    }
    // Saving sphere-to-planewall contact information
    for (const auto &pair1 : contactforce->getWallSphereContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "SphereToPlanewallContact"
                 << ", "
                 << info.particleId1 << ", "
                 << info.particleId2 << ", "
                 << info.normalForce.transpose() << ", "
                 << info.tangentialForce.transpose() << ", "
                 << info.tangentialDisplacement.transpose() << ", "
                 << info.previousTangentialVelocity.transpose() << "\n";
        }
    }
    // Saving sphere-to-fiber contact information
    for (const auto &pair1 : contactforce->getSphereFiberContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "SphereToFiberContact"
                 << ", "
                 << info.particleId1 << ", "
                 << info.particleId2 << ", "
                 << info.normalForce.transpose() << ", "
                 << info.tangentialForce.transpose() << ", "
                 << info.tangentialDisplacement.transpose() << ", "
                 << info.previousTangentialVelocity.transpose() << "\n";
        }
    }

    // Saving fiber-to-fiber contact information
    for (const auto &pair1 : contactforce->getFiberFiberContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "FiberToFiberContact"
                 << ", "
                 << info.particleId1 << ", "
                 << info.particleId2 << ", "
                 << info.normalForce.transpose() << ", "
                 << info.tangentialForce.transpose() << ", "
                 << info.tangentialDisplacement.transpose() << ", "
                 << info.previousTangentialVelocity.transpose() << "\n";
        }
    }

    // Saving wall-to-fiber contact information
    for (const auto &pair1 : contactforce->getWallFiberContactInformationList())
    {
        for (const auto &pair2 : pair1.second)
        {
            const ContactInformation &info = pair2.second;
            file << "WallToFiberContact"
                 << ", "
                 << info.particleId1 << ", " // Note: For wall contacts, one ID may be a wall ID or similar identifier.
                 << info.particleId2 << ", "
                 << info.normalForce.transpose() << ", "
                 << info.tangentialForce.transpose() << ", "
                 << info.tangentialDisplacement.transpose() << ", "
                 << info.previousTangentialVelocity.transpose() << "\n";
        }
    }

    file.close();
}
void DEMProperties::parseRandomParticle(std::istringstream &iss)
{

    std::string type;
    iss >> type;
    if (type == "SPHERE")
    {
        int categoryId, subTypeId, state, count;
        double xmin, xmax, ymin, ymax, zmin, zmax;
        iss >> categoryId >> subTypeId >> state >> count >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;

        int id_index = sphereparticles.size();

        gernerateSphereFlag = true;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> disX(xmin, xmax);
        std::uniform_real_distribution<> disY(ymin, ymax);
        std::uniform_real_distribution<> disZ(zmin, zmax);
        int failedSpheres = 0;
        double sphereRadius = manager->getSphereProperties(PropertyTypeID(categoryId, subTypeId))->getRadius();

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
                Eigen::Vector3d center = Eigen::Vector3d(x, y, z);
                for (const auto &plane : planewalls)
                {
                    // Retrieve the plane wall's point and normal
                    Eigen::Vector3d planePoint = plane->getCorner1();

                    Eigen::Vector3d planeNormal = plane->getNormal();

                    // Compute the vector from a point on the plane to the sphere's center
                    Eigen::Vector3d vecToSphereCenter = center - planePoint;

                    // Compute the distance from the sphere's center to the plane
                    double distanceToPlane = vecToSphereCenter.dot(planeNormal);

                    // Calculate the overlap (penetration depth)
                    double overlap_pw = sphereRadius - std::abs(distanceToPlane);
                    if (overlap_pw > 0)
                    {
                        validPosition = false;
                        break;
                    }
                }
                if (rectangularcontainer && validPosition)
                {
                    for (const auto &plane : rectangularcontainer->getPlaneWall())
                    {
                        // Retrieve the plane wall's point and normal
                        Eigen::Vector3d planePoint = plane->getCorner1();

                        Eigen::Vector3d planeNormal = plane->getNormal();

                        // Compute the vector from a point on the plane to the sphere's center
                        Eigen::Vector3d vecToSphereCenter = center - planePoint;

                        // Compute the distance from the sphere's center to the plane
                        double distanceToPlane = vecToSphereCenter.dot(planeNormal);

                        // Calculate the overlap (penetration depth)
                        double overlap_pw = sphereRadius - std::abs(distanceToPlane);
                        if (overlap_pw > 0)
                        {
                            validPosition = false;
                            break;
                        }
                    }
                }

                if (validPosition)
                {
                    for (const auto &existparticle : sphereparticles)
                    {
                        double distance = (existparticle->getPosition() - center).norm();
                        double sphere1Radius = manager->getSphereProperties(existparticle->getType())->getRadius();
                        double overlap_pp = (sphere1Radius + sphereRadius) - distance;

                        if (overlap_pp > 0)
                        {
                            validPosition = false;
                            break;
                        }
                    }
                }
                if (validPosition)
                {
                    for (const auto &bond : fiberbonds)
                    {
                        Eigen::Vector3d onenode1position = fibershpereparticles[bond->getNode1()]->getPosition();
                        Eigen::Vector3d onenode2position = fibershpereparticles[bond->getNode2()]->getPosition();

                        double fiberRadius = manager->getFiberProperties(bond->getType())->getRadius();
                        double t;
                        Eigen::Vector3d projection;
                        SphereCylinderBond::computeOverlap(onenode1position, onenode2position, center, t, projection);
                        double distance = (center - projection).norm();

                        double overlap = (sphereRadius + fiberRadius) - distance;
                        if (overlap > 0)
                        {
                            validPosition = false;
                            break;
                        }

                        if (!validPosition)
                            break;
                    }
                }
            } while (!validPosition && count_number < 100);
            if (validPosition)
            {
                auto sp = std::make_shared<SphereParticle>(id_index, PropertyTypeID(categoryId, subTypeId), state, manager, Eigen::Vector3d(x, y, z));
                sphereparticles.push_back(sp);
                id_index++;
            }
            else
            {
                gernerateSphereFlag = false;
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

                gernerateInfor[PropertyTypeID(categoryId, subTypeId)] = oss.str();

                break;
            }
        }
        if (gernerateSphereFlag)
        {
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

            gernerateInfor[PropertyTypeID(categoryId, subTypeId)] = oss.str();
        }
        else
        {
            std::cout << "Can't generate all spheres. Failed spheres: " << failedSpheres << std::endl;
        }
    }
    else if (type == "FIBER")
    {
        int categoryId, subTypeId, startstate, endstate, count;
        double xmin, xmax, ymin, ymax, zmin, zmax;
        iss >> categoryId >> subTypeId >> startstate >> endstate >> count >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;

        const auto &fiberproperties = manager->getFiberProperties(PropertyTypeID(categoryId, subTypeId));
        //FiberProperties copyfiberproperties = *fiberproperties;

        int bond_id_index = fiberbonds.size();
        int sphere_id_index = fibershpereparticles.size();
        int fiber_id_index = fibers.size();
        gernerateFiberFlag = true;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> disX(xmin, xmax);
        std::uniform_real_distribution<> disY(ymin, ymax);
        std::uniform_real_distribution<> disZ(zmin, zmax);

        std::uniform_real_distribution<> oriX(-1, 1);
        std::uniform_real_distribution<> oriY(-1, 1);
        std::uniform_real_distribution<> oriZ(-1, 1);
        int failedFibers = 0;
        double elementlength = fiberproperties->getElementlength();
        double fiberRadius = fiberproperties->getRadius();

        for (int i = 0; i < count; ++i)
        {
            double x, y, z, ox, oy, oz;
            bool validPosition = false;
            int count_number = 0;
            std::vector<Eigen::Vector3d> tempsphereposition;
            do
            {
                count_number++;
                validPosition = true;
                // Random starting position of the fiber
                x = disX(gen);
                y = disY(gen);
                z = disZ(gen);
                ox = oriX(gen);
                oy = oriY(gen);
                oz = oriZ(gen);

                Eigen::Vector3d startPos(x, y, z);
                // Random orientation of the fiber
                Eigen::Vector3d orientation(ox, oy, oz);
                orientation.normalize();
                for (int k = 0; k < fiberproperties->getNodenumber(); ++k)
                {

                    Eigen::Vector3d position = startPos + k * elementlength * orientation;

                    tempsphereposition.push_back(position);
                }

                for (const auto &plane : planewalls)
                {
                    // Retrieve the plane wall's point and normal
                    const Eigen::Vector3d &planePoint = plane->getCorner1();
                    const Eigen::Vector3d &planeNormal = plane->getNormal();
                    for (int j = 0; j < tempsphereposition.size() - 1; ++j)
                    {
                        Eigen::Vector3d projection;
                        SphereCylinderBond::computeOverlap(tempsphereposition[j], tempsphereposition[j + 1], plane, projection);
                        // Compute the vector from a point on the plane to the sphere's center
                        Eigen::Vector3d vecToSphereCenter = projection - planePoint;

                        // Compute the distance from the sphere's center to the plane
                        double distanceToPlane = vecToSphereCenter.dot(planeNormal);

                        // Calculate the overlap (penetration depth)
                        double overlap_pw = fiberRadius - std::abs(distanceToPlane);
                        if (overlap_pw > 0)
                        {
                            validPosition = false;
                            break;
                        }
                    }
                    if (!validPosition)
                        break;
                }
                if (rectangularcontainer && validPosition)
                {
                    for (const auto &plane : rectangularcontainer->getPlaneWall())
                    {
                        // Retrieve the plane wall's point and normal
                        const Eigen::Vector3d &planePoint = plane->getCorner1();
                        const Eigen::Vector3d &planeNormal = plane->getNormal();
                        for (int j = 0; j < tempsphereposition.size() - 1; ++j)
                        {
                            Eigen::Vector3d projection;
                            SphereCylinderBond::computeOverlap(tempsphereposition[j], tempsphereposition[j + 1], plane, projection);
                            // Compute the vector from a point on the plane to the sphere's center
                            Eigen::Vector3d vecToSphereCenter = projection - planePoint;

                            // Compute the distance from the sphere's center to the plane
                            double distanceToPlane = vecToSphereCenter.dot(planeNormal);

                            // Calculate the overlap (penetration depth)
                            double overlap_pw = fiberRadius - std::abs(distanceToPlane);
                            if (overlap_pw > 0)
                            {
                                validPosition = false;
                                break;
                            }
                        }
                        if (!validPosition)
                            break;
                    }
                }
                // Check overlap with existing spheres
                if (validPosition)
                {
                    for (const auto &sphere : sphereparticles)
                    {
                        // Retrieve the sphere's center and radius
                        Eigen::Vector3d sphereCenter = sphere->getPosition();
                        double sphereRadius = manager->getSphereProperties(sphere->getType())->getRadius();
                        for (int j = 0; j < tempsphereposition.size() - 1; ++j)
                        {
                            Eigen::Vector3d projection;
                            double t;
                            SphereCylinderBond::computeOverlap(tempsphereposition[j], tempsphereposition[j + 1], sphereCenter, t, projection);
                            double distance = (projection - sphereCenter).norm();

                            double overlap = (fiberRadius + sphereRadius) - distance;
                            if (overlap > 0)
                            {
                                validPosition = false;
                                break;
                            }
                        }
                        if (!validPosition)
                            break;
                    }
                }

                // Check overlap with existing fibers
                if (validPosition)
                {
                    for (const auto &bond : fiberbonds)
                    {
                        Eigen::Vector3d onenode1position = fibershpereparticles[bond->getNode1()]->getPosition();
                        Eigen::Vector3d onenode2position = fibershpereparticles[bond->getNode2()]->getPosition();

                        double fiber1Radius = manager->getFiberProperties(bond->getType())->getRadius();
                        for (int j = 0; j < tempsphereposition.size() - 1; ++j)
                        {
                            double s = 0, t = 0;
                            SphereCylinderBond::computeOverlap(tempsphereposition[j], tempsphereposition[j + 1], s,
                                                               onenode1position, onenode2position, t);
                            Eigen::Vector3d projection1 = (1 - s) * tempsphereposition[j] + s * tempsphereposition[j + 1];
                            Eigen::Vector3d projection2 = (1 - t) * onenode1position + t * onenode2position;
                            double distance = (projection1 - projection2).norm();

                            double overlap = (fiberRadius + fiber1Radius) - distance;
                            if (overlap > 0)
                            {
                                validPosition = false;
                                break;
                            }
                        }
                        if (!validPosition)
                            break;
                    }
                }

            } while (!validPosition && count_number < 100);
            if (validPosition)
            {
                auto fiber = std::make_shared<Fiber>(fiber_id_index, PropertyTypeID(categoryId, subTypeId), startstate, endstate, manager);
                std::vector<int> spherecylinderbondID;
                for (int k = 0; k < fiberproperties->getNodenumber(); ++k)
                {
                    auto createSphereParticle = [&](int state)
                    {
                        Eigen::Vector3d position = tempsphereposition[k];

                        return std::make_shared<SphereParticle>(
                            sphere_id_index + k, PropertyTypeID(categoryId, subTypeId), state, manager, position);
                    };

                    // Determine the state of the SphereParticle
                    int state = (k == 0) ? startstate : (k == fiberproperties->getNodenumber() - 1) ? endstate
                                                                                                       : 1;
                    auto sp = createSphereParticle(state);
                    fibershpereparticles.push_back(sp);

                    // Determine the SphereCylinderBond neighbors
                    int neighborElement1 = (k == 0) ? -1 : bond_id_index + k - 1;
                    int neighborElement2 = (k == fiberproperties->getNodenumber() - 2) ? -1 : bond_id_index + k + 1;
                    // Create SphereCylinderBond except for the last element
                    if (k != fiberproperties->getNodenumber() - 1)
                    {
                        auto bond = std::make_shared<SphereCylinderBond>(
                            bond_id_index + k,
                            PropertyTypeID(categoryId, subTypeId),
                            manager,
                            fiber_id_index,
                            sphere_id_index + k,
                            sphere_id_index + k + 1,
                            neighborElement1,
                            neighborElement2);
                        fiberbonds.push_back(bond);
                        spherecylinderbondID.push_back(bond_id_index + k);
                    }
                }
                fiber->setSphereCylinderBond(spherecylinderbondID);
                fibers.push_back(fiber);
                sphere_id_index += fiberproperties->getNodenumber();
                fiber_id_index++;
                bond_id_index += fiberproperties->getNodenumber() - 1;
            }
            else
            {
                gernerateFiberFlag = false;
                failedFibers = count - i;
                // Use an ostringstream to format the data as a string
                std::ostringstream oss;
                oss << type << " "
                    << categoryId << " "
                    << subTypeId << " "
                    << startstate << " "
                    << endstate << " "
                    << failedFibers << " "
                    << xmin << " "
                    << xmax << " "
                    << ymin << " "
                    << ymax << " "
                    << zmin << " "
                    << zmax;

                gernerateInfor[PropertyTypeID(categoryId, subTypeId)] = oss.str();

                break;
            }
        }
        if (gernerateFiberFlag)
        {
            // Use an ostringstream to format the data as a string
            std::ostringstream oss;
            oss << type << " "
                << categoryId << " "
                << subTypeId << " "
                << startstate << " "
                << endstate << " "
                << failedFibers << " "
                << xmin << " "
                << xmax << " "
                << ymin << " "
                << ymax << " "
                << zmin << " "
                << zmax;

            gernerateInfor[PropertyTypeID(categoryId, subTypeId)] = oss.str();
        }
        else
        {
            std::cout << "Can't generate all fibers. Failed fiber: " << failedFibers << std::endl;
        }
    }
}
void DEMProperties::parseSpecificParticle(std::istringstream &iss)
{
    std::string type;
    iss >> type;
    if (type == "SPHERE")
    {
        int id, categoryId, subTypeId, state;
        double x, y, z, vx, vy, vz;
        iss >> id >> categoryId >> subTypeId >> state >> x >> y >> z >> vx >> vy >> vz;
        // Create and add the specific particle to your simulation
        auto sp = std::make_shared<SphereParticle>(id, PropertyTypeID(categoryId, subTypeId), state, manager, Eigen::Vector3d(x, y, z), Eigen::Vector3d(vx, vy, vz));

        sphereparticles.push_back(sp);
    }
    else if (type == "FIBER")
    {
        int id, categoryId, subTypeId, startstate, endstate;
        double startx, starty, startz, vx, vy, vz;
        iss >> id >> categoryId >> subTypeId >> startstate >> endstate >> startx >> starty >> startz >> vx >> vy >> vz;
        std::string lauout;
        iss >> lauout;
        auto fiber = std::make_shared<Fiber>(id, PropertyTypeID(categoryId, subTypeId), startstate, endstate, manager);
        const auto &fiberproperties = manager->getFiberProperties(PropertyTypeID(categoryId, subTypeId));
        FiberProperties copyfiberproperties = *fiberproperties;
        int bond_id_index = fiberbonds.size();
        int sphere_id_index = fibershpereparticles.size();
        std::vector<int> spherecylinderbondID;
        for (int k = 0; k < copyfiberproperties.getNodenumber(); ++k)
        {
            // Lambda function to create a SphereParticle based on layout
            auto createSphereParticle = [&](int state, const std::string &layout)
            {
                Eigen::Vector3d position;
                if (layout == "X")
                {
                    position = Eigen::Vector3d(startx + k * copyfiberproperties.getElementlength(), starty, startz);
                }
                else if (layout == "Y")
                {
                    position = Eigen::Vector3d(startx, starty + k * copyfiberproperties.getElementlength(), startz);
                }
                else if (layout == "Z")
                {
                    position = Eigen::Vector3d(startx, starty, startz + k * copyfiberproperties.getElementlength());
                }
                return std::make_shared<SphereParticle>(
                    sphere_id_index + k, PropertyTypeID(categoryId, subTypeId), state, manager, position, Eigen::Vector3d(vx, vy, vz));
            };

            // Determine the state of the SphereParticle
            int state = (k == 0) ? startstate : (k == copyfiberproperties.getNodenumber() - 1) ? endstate
                                                                                               : 1;
            auto sp = createSphereParticle(state, lauout);
            fibershpereparticles.push_back(sp);

            // Determine the SphereCylinderBond neighbors
            int neighborElement1 = (k == 0) ? -1 : bond_id_index + k - 1;
            int neighborElement2 = (k == copyfiberproperties.getNodenumber() - 2) ? -1 : bond_id_index + k + 1;

            // Create SphereCylinderBond except for the last element
            if (k != copyfiberproperties.getNodenumber() - 1)
            {
                auto bond = std::make_shared<SphereCylinderBond>(
                    bond_id_index + k,
                    PropertyTypeID(categoryId, subTypeId),
                    manager,
                    id,
                    sphere_id_index + k,
                    sphere_id_index + k + 1,
                    neighborElement1,
                    neighborElement2);
                fiberbonds.push_back(bond);
                spherecylinderbondID.push_back(bond_id_index + k);
            }
        }
        fiber->setSphereCylinderBond(spherecylinderbondID);
        fibers.push_back(fiber);
    }
}
void DEMProperties::generateRemainingParticles()
{
    if (!gernerateFiberFlag || !gernerateSphereFlag)
    {

        for (const auto &sub : gernerateInfor)
        {

            std::istringstream iss(sub.second);
            parseRandomParticle(iss);
        }
    }
}
double DEMProperties::getAverageVelocity()
{
    double totalVelocity = 0;

    for (const auto &sphere : sphereparticles)
    {
        totalVelocity += sphere->getVelocity().norm();
    }
    return totalVelocity / sphereparticles.size();
}
void DEMProperties::initialSimulation()
{
    contactforce->addParticleProperties(manager);
    double maxRadius = std::numeric_limits<double>::min();
    double maxYoungmodulus = std::numeric_limits<double>::min();
    double maxBondmodulus = std::numeric_limits<double>::min();
    double minBondlength = std::numeric_limits<double>::max();
    double mindensity = std::numeric_limits<double>::max();
    double minRadius = std::numeric_limits<double>::max();
    double mindensityL = std::numeric_limits<double>::max();


    for (const auto &particleproperties : manager->getParticleProperties())
    {
        if (auto& sphereProps = std::dynamic_pointer_cast<SphereProperties>(particleproperties.second))
        {
            double radius = sphereProps->getRadius();
            if (radius > maxRadius)
            {
                maxRadius = radius;
            }
            if(radius < minRadius)
            {
                minRadius = radius;
            }
            double Youngmodulus = sphereProps->getYoungModulus();
            if (Youngmodulus > maxYoungmodulus)
            {
                maxYoungmodulus = Youngmodulus;
            }
            double density = sphereProps->getDensity();
            if (density < mindensity)
            {
                mindensity = density;
            }
        }
        else if (auto& fiberProps = std::dynamic_pointer_cast<FiberProperties>(particleproperties.second))
        {
            double radius = fiberProps->getRadius();
            if (radius > maxRadius)
            {
                maxRadius = radius;
            }
            if(radius < minRadius)
            {
                minRadius = radius;
            }
            double Youngmodulus = fiberProps->getYoungModulus();
            if (Youngmodulus > maxYoungmodulus)
            {
                maxYoungmodulus = Youngmodulus;
            }
            double Bondmodulus = fiberProps->getNormalmodulus();
            if(Bondmodulus > maxBondmodulus)
            {
                maxBondmodulus = Bondmodulus;
            }
            Bondmodulus = fiberProps->getShearmodulus();
            if(Bondmodulus > maxBondmodulus)
            {
                maxBondmodulus = Bondmodulus;
            }
            Bondmodulus = fiberProps->getTwistmodulus();
            if(Bondmodulus > maxBondmodulus)
            {
                maxBondmodulus = Bondmodulus;
            }
            Bondmodulus = fiberProps->getBendingmodulus();
            if(Bondmodulus > maxBondmodulus)
            {
                maxBondmodulus = Bondmodulus;
            }


            double Bondlength = fiberProps->getElementlength();
            if(Bondlength < minBondlength)
            {
                minBondlength = Bondlength;
            }

            double densityL = fiberProps->getNodemass()  / Bondlength;
            if(densityL < mindensityL)
            {
                mindensityL = densityL;
            }

            double density = fiberProps->getDensity();
            if (density < mindensity)
            {
                mindensity = density;
            }
        }
    }


    


   

    double temp_t = 0.1 * minBondlength * sqrt(mindensity / (maxBondmodulus*PI*maxRadius*maxRadius));
    if (timestep > temp_t)
    {
        timestep = temp_t;
    }
    temp_t = (0.1 * 2*PI*minRadius) / (0.1631 * 0.3 + 0.8766) * sqrt(1.3 * mindensity / (2*maxYoungmodulus));
    if (timestep > temp_t)
    {
        timestep = temp_t;
    }
    gridbasedcontactdetection->initial(simulationdimensions.x(), simulationdimensions.y(), simulationdimensions.z(), minRadius);
}
void DEMProperties::handleCollisions()
{
    std::unordered_map<int, std::unordered_set<int>> sphere_sphere_contact_paris;
    std::unordered_map<int, std::unordered_set<int>> sphere_fiber_contact_paris;

    std::unordered_map<int, std::unordered_set<int>> planewall_sphere_contact_paris;
    std::unordered_map<int, std::unordered_set<int>> planewall_fiber_contact_paris;

    std::unordered_map<int, std::unordered_set<int>> fiber_fiber_contact_paris;

    std::unordered_map<int, std::unordered_set<int>> rectangularcontainer_sphere_contact_paris;
    std::unordered_map<int, std::unordered_set<int>> rectangularcontainer_fiber_contact_paris;

    gridbasedcontactdetection->ParticleBroadPhase(sphereparticles, fiberbonds, fibershpereparticles,
                                                  sphere_sphere_contact_paris, sphere_fiber_contact_paris, fiber_fiber_contact_paris);
    for (const auto &contact_list : sphere_sphere_contact_paris)
    {
        auto id1 = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            auto &sp1 = sphereparticles[id1];
            auto &sp2 = sphereparticles[id2];

            contactforce->computeSphereSphereForce(sp1, sp2, timestep);
        }
    }

    for (const auto &contact_list : sphere_fiber_contact_paris)
    {
        auto id1 = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            auto &sp1 = sphereparticles[id1];
            auto &sc2 = fiberbonds[id2];

            auto node1 = fibershpereparticles[sc2->getNode1()];
            auto node2 = fibershpereparticles[sc2->getNode2()];

            contactforce->computeSphereFiberForce(sp1, sc2, node1, node2, timestep);
        }
    }

    for (const auto &contact_list : fiber_fiber_contact_paris)
    {
        auto id1 = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            const auto &sc1 = fiberbonds[id1];
            const auto &sc2 = fiberbonds[id2];

            auto &sconenode1 = fibershpereparticles[sc1->getNode1()];
            auto &sconenode2 = fibershpereparticles[sc1->getNode2()];
            auto &sctwonode1 = fibershpereparticles[sc2->getNode1()];
            auto &sctwonode2 = fibershpereparticles[sc2->getNode2()];

            contactforce->computeFiberFiberForce(sc1, sconenode1, sconenode2, sc2, sctwonode1, sctwonode2, timestep);
        }
    }

    gridbasedcontactdetection->planewallBroadPhase(planewalls, planewall_sphere_contact_paris, planewall_fiber_contact_paris);
    for (const auto &contact_list : planewall_sphere_contact_paris)
    {
        auto wall_id = contact_list.first;
        for (auto id : contact_list.second)
        {
            auto &wall = planewalls[wall_id];
            auto &sphere = sphereparticles[id];

            contactforce->computePlaneWallSphereForce(wall, sphere, timestep);
        }
    }
    for (const auto &contact_list : planewall_fiber_contact_paris)
    {
        auto wall_id = contact_list.first;
        for (auto id2 : contact_list.second)
        {
            auto &wall = planewalls[wall_id];
            auto &sc2 = fiberbonds[id2];

            auto node1 = fibershpereparticles[sc2->getNode1()];
            auto node2 = fibershpereparticles[sc2->getNode2()];

            contactforce->computePlaneWallFiberForce(wall, sc2, node1, node2, timestep);
        }
    }

    if (rectangularcontainer)
    {
        const auto &rcplanewalls = rectangularcontainer->getPlaneWall();
        gridbasedcontactdetection->RectangularContainerBroadPhase(rectangularcontainer, rectangularcontainer_sphere_contact_paris, rectangularcontainer_fiber_contact_paris);
        for (const auto &contact_list : rectangularcontainer_sphere_contact_paris)
        {
            auto wall_id = contact_list.first;
            for (auto id : contact_list.second)
            {
                auto &wall = rcplanewalls[wall_id];
                auto &sphere = sphereparticles[id];

                contactforce->computePlaneWallSphereForce(wall, sphere, timestep);
            }
        }
        for (const auto &contact_list : rectangularcontainer_fiber_contact_paris)
        {
            auto wall_id = contact_list.first;
            for (auto id2 : contact_list.second)
            {
                auto &wall = rcplanewalls[wall_id];
                auto &sc2 = fiberbonds[id2];

                auto node1 = fibershpereparticles[sc2->getNode1()];
                auto node2 = fibershpereparticles[sc2->getNode2()];

                contactforce->computePlaneWallFiberForce(wall, sc2, node1, node2, timestep);
            }
        }
    }
}
void DEMProperties::bondforce()
{
     for (const auto &particle : fiberbonds)
    {
        auto& node1 = fibershpereparticles[particle->getNode1()];
        auto& node2 = fibershpereparticles[particle->getNode2()];
        particle->updateBond(node1,node2,timestep);

    }
}
void DEMProperties::applyExternalForces()
{
    for (auto &sphere : sphereparticles)
    {

        sphere->addForce(globalforces);
    }
    for (auto &indexforce : specificforces)
    {
        sphereparticles[indexforce.first]->addForce(indexforce.second);
    }
    // add gravity force
}
void DEMProperties::motion()
{
    // motion particles
    for (auto &sphere : sphereparticles)
    {
        if (sphere->getState() == 1)
        {
            sphere->updateVelocity(timestep, gravity);
            sphere->updateOmega(timestep);
            sphere->updatePosition(timestep);
        }
        sphere->resetForce();
        sphere->resetTorque();
    }

    for (auto &sphere : fibershpereparticles)
    {
        if (sphere->getState() == 1)
        {
            sphere->updateVelocity(timestep, gravity);
            sphere->updateOmega(timestep);
            sphere->updatePosition(timestep);
        }
        sphere->resetForce();
        sphere->resetTorque();
    }
    
    
    
    // motion rectangularcontainer
    if (rectangularcontainer)
    {

        Eigen::Vector3d resultantForce = Eigen::Vector3d::Zero();
        for (auto &planewall : rectangularcontainer->getPlaneWall())
        {
            resultantForce += planewall->getForce();
            planewall->resetForce();
        }
        rectangularcontainer->addForce(resultantForce);

        if (rectangularcontainer->getState() == 1)
        {
            double mass = rectangularcontainer->getMass();
            Eigen::Vector3d velocity = rectangularcontainer->getVelocity();
            velocity += (rectangularcontainer->getForce() / mass + gravity) * timestep;
            rectangularcontainer->setVelocity(velocity);
            Eigen::Vector3d displancement = velocity * timestep;
            Eigen::Vector3d newCenter = rectangularcontainer->getCenter() + displancement;
            rectangularcontainer->setCenter(newCenter);

            for (auto &planewall : rectangularcontainer->getPlaneWall())
            {
                Eigen::Vector3d p1 = planewall->getCorner1();
                Eigen::Vector3d p2 = planewall->getCorner2();
                Eigen::Vector3d p3 = planewall->getCorner3();
                p1 += displancement;
                p2 += displancement;
                p3 += displancement;
                planewall->setCorner1(p1);
                planewall->setCorner2(p2);
                planewall->setCorner3(p3);
                planewall->setVelocity(velocity);
            }
        }

        rectangularcontainer->resetForce();
    }


    contactforce->updateContactInformation();
}
void DEMProperties::initial_task()
{
    for (auto &subtask : dynamicBoundaryTask)
    {
        if (subtask.first == "RECTANGULARCONTAINER")
        {

            if (auto &springoscillator = std::dynamic_pointer_cast<SpringOscillator>(subtask.second[0]))
            {
                double y0 = rectangularcontainer->getCenter().y();
                double position = springoscillator->getequilibriumPosition();
                position += y0;
                springoscillator->setequilibriumPosition(position);
                rectangularcontainer->setMass(manager);

                rectangularcontainer->setState(1);
                for (auto &planewall : rectangularcontainer->getPlaneWall())
                {
                    planewall->setState(1);
                }
            }
        }
    }
}
void DEMProperties::dotask()
{
    for (auto &subtask : dynamicBoundaryTask)
    {
        if (subtask.first == "RECTANGULARCONTAINER")
        {
            Eigen::Vector3d resultantForce = Eigen::Vector3d::Zero();
            Eigen::Vector3d v = rectangularcontainer->getVelocity();
            if (auto &springoscillator = std::dynamic_pointer_cast<SpringOscillator>(subtask.second[0]))
            {
                resultantForce = springoscillator->getForce(rectangularcontainer->getCenter().y(), rectangularcontainer->getVelocity().y());
            }
            rectangularcontainer->addForce(resultantForce);
        }
    }
}
