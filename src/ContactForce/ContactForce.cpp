#include "ContactForce.h"

void ContactForce::addParticleProperties(const std::shared_ptr<ParticlePropertyManager> manager)
{
    for (auto &sub_propertys : manager->getParticleProperties())
    {
        PropertyTypeID i = sub_propertys.first;
        const double radius_i = sub_propertys.second->getRadius();
        const double youngs_modulus_i = sub_propertys.second->getYoungModulus();
        const double poisson_ratio_i = sub_propertys.second->getPoissonRatio();
        const double restitution_i = sub_propertys.second->getRestitution();
        const double contact_damping_coefficient_i = -log(restitution_i) / sqrt(PI * PI + log(restitution_i) * log(restitution_i));
        const double sliding_friction_i = sub_propertys.second->getSlidingFriction();
        const double rolling_friction_i = sub_propertys.second->getRollingFriction();
        const double mass_i = sub_propertys.second->getMass();

        for (auto &another_sub_propertys : manager->getParticleProperties())
        {
            PropertyTypeID j = another_sub_propertys.first;
            const double radius_j = another_sub_propertys.second->getRadius();
            const double youngs_modulus_j = another_sub_propertys.second->getYoungModulus();
            const double poisson_ratio_j = another_sub_propertys.second->getPoissonRatio();
            const double restitution_j = another_sub_propertys.second->getRestitution();
            const double contact_damping_coefficient_j = -log(restitution_j) / sqrt(PI * PI + log(restitution_j) * log(restitution_j));
            const double sliding_friction_j = another_sub_propertys.second->getSlidingFriction();
            const double rolling_friction_j = another_sub_propertys.second->getRollingFriction();
            const double mass_j = another_sub_propertys.second->getMass();

            effectiveradius[i][j] = 2 * (radius_i * radius_j) / (radius_i + radius_j);

            effectivemass[i][j] = 2 * (mass_i * mass_j) / (mass_i + mass_j);

            effectiveyoungsmodulus[i][j] = (youngs_modulus_i * youngs_modulus_j) /
                                           ((youngs_modulus_j * (1 - poisson_ratio_i * poisson_ratio_i)) +
                                            (youngs_modulus_i * (1 - poisson_ratio_j * poisson_ratio_j)));

            effectiveshearmodulus[i][j] = (youngs_modulus_i * youngs_modulus_j) /
                                          (2 * ((youngs_modulus_j * (2 - poisson_ratio_i) *
                                                 (1 + poisson_ratio_i)) +
                                                (youngs_modulus_i * (2 - poisson_ratio_j) *
                                                 (1 + poisson_ratio_j))));

            modelparameterbeta[i][j] = 2 * (contact_damping_coefficient_i * contact_damping_coefficient_j) /
                                       (contact_damping_coefficient_i + contact_damping_coefficient_j + DBL_MIN);

            effectiveslidingfriction[i][j] = 2 * sliding_friction_i * sliding_friction_j /
                                             (sliding_friction_i + sliding_friction_j + DBL_MIN);

            effectiverollingfriction[i][j] = 2 * rolling_friction_i * rolling_friction_j /
                                             (rolling_friction_i + rolling_friction_j + DBL_MIN);
        }
    }
}

void ContactForce::computeSphereSphereForce(std::shared_ptr<SphereParticle> sphere1, std::shared_ptr<SphereParticle> sphere2, double timeStep)
{

    auto manager = sphere1->getParticlePropertyManager();

    PropertyTypeID particletype_one = sphere1->getType();
    PropertyTypeID particletype_two = sphere2->getType();

    auto radius1 = manager->getSphereProperties(particletype_one)->getRadius();
    auto radius2 = manager->getSphereProperties(particletype_two)->getRadius();

    // Compute the distance between the centers of the two spheres
    Eigen::Vector3d position1 = sphere1->getPosition();
    Eigen::Vector3d position2 = sphere2->getPosition();
    double distance = (position1 - position2).norm();

    // Compute the overlap (positive if spheres intersect, zero or negative otherwise)
    double normal_overlap = (radius1 + radius2) - distance;

    // The Hertz-Mindlin  contact model

    if (normal_overlap > 0)
    {
        Eigen::Vector3d normal_unit_vector = position2 - position1;
        normal_unit_vector.normalize();

        Eigen::Vector3d v1 = sphere1->getVelocity();

        Eigen::Vector3d v2 = sphere2->getVelocity();

        Eigen::Vector3d w1 = sphere1->getOmega();

        Eigen::Vector3d w2 = sphere2->getOmega();

        Eigen::Vector3d contact_relative_velocity = (v1 - v2) + (radius1 * w1 + radius2 * w2).cross(normal_unit_vector);

        double normal_relative_velocity_value = contact_relative_velocity.dot(normal_unit_vector);

        Eigen::Vector3d normal_v = normal_relative_velocity_value * normal_unit_vector;

        // Calculation of normal relative velocity
        Eigen::Vector3d normal_relative_velocity = normal_relative_velocity_value * normal_unit_vector;
        // Calculation of tangential relative velocity
        Eigen::Vector3d tangential_relative_velocity = contact_relative_velocity - normal_relative_velocity;
        // modify tangential_overlap
        int id1 = sphere1->getId();
        int id2 = sphere2->getId();

        auto &innerMap = contactinformationlist[id1]; // Access or create the inner map
        auto innerIt = innerMap.find(id2);

        if (innerIt == innerMap.end())
        {
            // Contact information for this pair does not exist, create a new one
            ContactInformation newContactInfo(id1, id2);
            innerMap[id2] = newContactInfo; // Insert the new contact information
        }

        ContactInformation &contact_info = contactinformationlist[id1][id2];

        Eigen::Vector3d modified_tangential_displacement = contact_info.tangentialDisplacement +  contact_info.previousTangentialVelocity * timeStep;
        // Updating the contact_info container based on the new calculated values
        contact_info.tangentialDisplacement = modified_tangential_displacement;
        contact_info.previousTangentialVelocity = tangential_relative_velocity;

        double radius_times_overlap_sqrt = sqrt(effectiveradius[particletype_one][particletype_two] * normal_overlap);

        double sn = 2 * effectiveyoungsmodulus[particletype_one][particletype_two] * radius_times_overlap_sqrt;

        double normal_contact_stiffness = 0.666666 * sn;

        double normal_damping_constant = 1.8257 * modelparameterbeta[particletype_one][particletype_two] * sqrt(sn * effectivemass[particletype_one][particletype_two]);

        // tangential stiffness
        double tangential_contact_stiffness = 8 * effectiveshearmodulus[particletype_one][particletype_two] * radius_times_overlap_sqrt;

        double tangential_damping_constant = 1.8257 * modelparameterbeta[particletype_one][particletype_two] * sqrt(tangential_contact_stiffness * effectivemass[particletype_one][particletype_two]);

        // Calculation of normal force using spring and dashpot normal forces
        Eigen::Vector3d normal_force =
            ((normal_contact_stiffness * normal_overlap) * normal_unit_vector) +
            ((normal_damping_constant * normal_relative_velocity_value) *
             normal_unit_vector);

        // Calculation of tangential force using spring and dashpot tangential
        // forces. Since we need dashpot tangential force in the gross sliding again,
        // we define it as a separate variable
        Eigen::Vector3d dashpot_tangential_force = tangential_damping_constant * tangential_relative_velocity;
        Eigen::Vector3d tangential_force = (tangential_contact_stiffness * modified_tangential_displacement) + dashpot_tangential_force;

        double coulomb_threshold = effectiveslidingfriction[particletype_one][particletype_one] * normal_force.norm();

        double tangential_force_value = tangential_force.norm();
        // Check for gross sliding
        if (tangential_force_value > coulomb_threshold)
        {
            // Gross sliding occurs and the tangential overlap and tangnetial
            // force are limited to Coulumb's criterion
            Eigen::Vector3d tangential_overlap = (coulomb_threshold * (tangential_force / (tangential_force_value + DBL_MIN)) - dashpot_tangential_force) /
                                                 (tangential_contact_stiffness + DBL_MIN);

            contact_info.tangentialDisplacement = tangential_overlap;
            tangential_force = (tangential_contact_stiffness * tangential_overlap) + dashpot_tangential_force;
        }

        // Calculation of torque
        // Torque caused by tangential force (tangential_torque)
        Eigen::Vector3d tangential_torque1 = (radius1 * normal_unit_vector).cross(tangential_force);
        Eigen::Vector3d tangential_torque2 = (radius2 * normal_unit_vector).cross(tangential_force);

        Eigen::Vector3d results_force = normal_force + tangential_force;

        sphere1->addForce(-results_force);
        sphere2->addForce(results_force);

        sphere1->addTorque(-tangential_torque1);
        sphere2->addTorque(-tangential_torque2);
        contact_info.normalForce = -normal_force;
        contact_info.tangentialForce = -tangential_force;
        nextcontactinformationlist[id1][id2] = contact_info;

    }
}

void ContactForce::computePlaneWallSphereForce(std::shared_ptr<PlaneWall> planewall, std::shared_ptr<SphereParticle> sphere, double timeStep)
{
    // Retrieve the sphere's center and radius
    Eigen::Vector3d sphereCenter = sphere->getPosition();

    auto manager = sphere->getParticlePropertyManager();

    PropertyTypeID spheretype = sphere->getType();
    PropertyTypeID walltype = planewall->getType();

    double sphereRadius = manager->getSphereProperties(spheretype)->getRadius();

    // Retrieve the plane wall's point and normal
    Eigen::Vector3d planePoint = planewall->getCorner1();

    Eigen::Vector3d planeNormal = planewall->getNormal();

    // Compute the vector from a point on the plane to the sphere's center
    Eigen::Vector3d vecToSphereCenter = sphereCenter - planePoint;

    // Compute the distance from the sphere's center to the plane
    double distanceToPlane = vecToSphereCenter.dot(planeNormal);

    if (distanceToPlane < 0)
    {
        planeNormal = -planeNormal;
        distanceToPlane = -distanceToPlane;
    }

    double normal_overlap = sphereRadius - distanceToPlane;

    if (normal_overlap > 0)
    {
        // Using contact_vector, the contact normal vector is obtained
        Eigen::Vector3d normal_unit_vector = -planeNormal;
        // Defining velocities and angular velocities of sphere
        Eigen::Vector3d sphere_velocity = sphere->getVelocity();

        Eigen::Vector3d sphere_omega = sphere->getOmega();
        // Calculation of contact relative velocity
        Eigen::Vector3d contact_relative_velocity = sphere_velocity + (sphereRadius * sphere_omega).cross(normal_unit_vector);

        double normal_relative_velocity_value = contact_relative_velocity.dot(normal_unit_vector);
        // Calculation of normal relative velocity
        Eigen::Vector3d normal_relative_velocity = normal_relative_velocity_value * normal_unit_vector;
        // Calculation of tangential relative velocity
        Eigen::Vector3d tangential_relative_velocity = contact_relative_velocity - normal_relative_velocity;
        // modify tangential_overlap
        int wallId = planewall->getId();
        int sphereId = sphere->getId();

        auto &innerMap = contactinformationlist[wallId]; // Access or create the inner map
        auto innerIt = innerMap.find(sphereId);

        if (innerIt == innerMap.end())
        {
            // Contact information for this pair does not exist, create a new one
            ContactInformation newContactInfo(wallId, sphereId);
            innerMap[sphereId] = newContactInfo; // Insert the new contact information
        }

        ContactInformation &contact_info = contactinformationlist[wallId][sphereId];

        Eigen::Vector3d modified_tangential_displacement = contact_info.tangentialDisplacement + contact_info.previousTangentialVelocity * timeStep;
        // Updating the contact_info container based on the new calculated values
        contact_info.tangentialDisplacement = modified_tangential_displacement;
        contact_info.previousTangentialVelocity = tangential_relative_velocity;

        double radius_times_overlap_sqrt = sqrt(manager->getSphereProperties(spheretype)->getRadius() * normal_overlap);

        double sn = 2 * effectiveyoungsmodulus[walltype][spheretype] * radius_times_overlap_sqrt;

        double normal_contact_stiffness = 0.666666 * sn;

        double normal_damping_constant = 1.8257 * modelparameterbeta[walltype][spheretype] * sqrt(sn * manager->getSphereProperties(spheretype)->getMass());

        // tangential stiffness
        double tangential_contact_stiffness = 8 * effectiveshearmodulus[walltype][spheretype] * radius_times_overlap_sqrt;

        double tangential_damping_constant = 1.8257 * modelparameterbeta[walltype][spheretype] * sqrt(tangential_contact_stiffness * manager->getSphereProperties(spheretype)->getMass());

        // Calculation of normal force using spring and dashpot normal forces
        Eigen::Vector3d normal_force =
            ((normal_contact_stiffness * normal_overlap) * normal_unit_vector) +
            ((normal_damping_constant * normal_relative_velocity_value) *
             normal_unit_vector);

        // Calculation of tangential force using spring and dashpot tangential
        // forces. Since we need dashpot tangential force in the gross sliding again,
        // we define it as a separate variable
        Eigen::Vector3d dashpot_tangential_force = tangential_damping_constant * tangential_relative_velocity;
        Eigen::Vector3d tangential_force = (tangential_contact_stiffness * modified_tangential_displacement) + dashpot_tangential_force;

        double coulomb_threshold = effectiveslidingfriction[walltype][spheretype] * normal_force.norm();

        double tangential_force_value = tangential_force.norm();
        // Check for gross sliding
        if (tangential_force_value > coulomb_threshold)
        {
            // Gross sliding occurs and the tangential overlap and tangnetial
            // force are limited to Coulumb's criterion
            Eigen::Vector3d tangential_overlap = (coulomb_threshold * (tangential_force / (tangential_force_value + DBL_MIN)) - dashpot_tangential_force) /
                                                 (tangential_contact_stiffness + DBL_MIN);

            contact_info.tangentialDisplacement = tangential_overlap;
            tangential_force = (tangential_contact_stiffness * tangential_overlap) + dashpot_tangential_force;
        }

        // Calculation of torque
        // Torque caused by tangential force (tangential_torque)
        Eigen::Vector3d tangential_torque1 = (sphereRadius * normal_unit_vector).cross(tangential_force);

        Eigen::Vector3d results_force = normal_force + tangential_force;

        sphere->addForce(-results_force);
        sphere->addTorque(-tangential_torque1);
        planewall->addForce(results_force);

        contact_info.normalForce = normal_force;
        contact_info.tangentialForce = tangential_force;
        nextcontactinformationlist[wallId][sphereId] = contact_info;
    }
}

void ContactForce::updateContactInformation()
{
    contactinformationlist = nextcontactinformationlist;
    nextcontactinformationlist.clear();

    wallcontactinformationlist = nextwallcontactinformationlist;
    nextwallcontactinformationlist.clear();
}