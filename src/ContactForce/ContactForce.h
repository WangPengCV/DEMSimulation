#pragma once
#include <memory>
#include <Eigen/Dense>
#include "SphereParticle.h"
<<<<<<< HEAD
#include "SphereCylinderBond.h"
#include "PlanewallProperties.h"
=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
#include "PlaneWall.h"
#include "ParticlePropertyManager.h"
#include <unordered_map>

class ContactInformation
{
public:
    int particleId1;                            // ID of the first particle in contact
    int particleId2;                            // ID of the second particle in contact
    Eigen::Vector3d tangentialDisplacement;     // Cumulative tangential displacement
    Eigen::Vector3d previousTangentialVelocity; // Tangential velocity in the previous timestep
    Eigen::Vector3d normalForce;                // Computed normal force during the last interaction
    Eigen::Vector3d tangentialForce;            // Computed tangential force during the last interaction

    ContactInformation(int id1, int id2) : particleId1(id1), particleId2(id2),
                                           tangentialDisplacement(Eigen::Vector3d::Zero()),
                                           previousTangentialVelocity(Eigen::Vector3d::Zero()),
                                           normalForce(Eigen::Vector3d::Zero()),
                                           tangentialForce(Eigen::Vector3d::Zero()) {}
<<<<<<< HEAD
    ContactInformation() : particleId1(-1), particleId2(-1), // -1 or another invalid ID
                           tangentialDisplacement(Eigen::Vector3d::Zero()),
                           previousTangentialVelocity(Eigen::Vector3d::Zero()),
                           normalForce(Eigen::Vector3d::Zero()),
                           tangentialForce(Eigen::Vector3d::Zero())
    {
    }
=======
    ContactInformation():particleId1(-1), particleId2(-1),  // -1 or another invalid ID
      tangentialDisplacement(Eigen::Vector3d::Zero()),
      previousTangentialVelocity(Eigen::Vector3d::Zero()),
      normalForce(Eigen::Vector3d::Zero()),
      tangentialForce(Eigen::Vector3d::Zero()) {}
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
};

class ContactForce
{
public:
<<<<<<< HEAD
    void addParticleProperties(const std::shared_ptr<ParticlePropertyManager> &manager);

    void computerSphereSphereEffective(std::shared_ptr<SphereProperties> &sphereproperties1, PropertyTypeID &type1,
                                       std::shared_ptr<SphereProperties> &sphereproperties2, PropertyTypeID &type2);

    void computerSphereFiberEffective();

    void computerSpherePlanewallEffective(std::shared_ptr<PlanewallProperties> &planewallproperties, PropertyTypeID &type1,
                                          std::shared_ptr<SphereProperties> &sphereproperties, PropertyTypeID &type2);

    void computerFiberPlanewallEffective();

    // Method to compute force between two spherical particles
    void computeSphereSphereForce(const std::shared_ptr<SphereParticle> &sphere1, const std::shared_ptr<SphereParticle> &sphere2, double timeStep);
    // Method to compute force between a sphere particle and a plane wall
    void computePlaneWallSphereForce(const std::shared_ptr<PlaneWall> &planeWall, const std::shared_ptr<SphereParticle> &sphere, double timeStep);
    // Method to compute force between two spherical particles
    void computeSphereFiberForce(const std::shared_ptr<SphereParticle> &sphere, const std::shared_ptr<SphereCylinderBond> & fiberbond,
                                 const std::shared_ptr<SphereParticle> &node1, const std::shared_ptr<SphereParticle> &node2, double timeStep);
    // Method to compute force between a sphere particle and a plane wall
    void computePlaneWallSphereForce(const std::shared_ptr<PlaneWall> &planeWall, const std::shared_ptr<SphereParticle> &sphere, double timeStep);
    // Method to compute force between a sphere particle and a plane wall
    void computePlaneWallSphereForce(const std::shared_ptr<PlaneWall> &planeWall, const std::shared_ptr<SphereParticle> &sphere, double timeStep);

    void updateContactInformation();

    std::unordered_map<int, std::unordered_map<int, ContactInformation>> &getContactInformationList();
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> &getWallContactInformationList();

private:
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> contactinformationlist;
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> nextcontactinformationlist;

    std::unordered_map<int, std::unordered_map<int, ContactInformation>> wallcontactinformationlist;
    std::unordered_map<int, std::unordered_map<int, ContactInformation>> nextwallcontactinformationlist;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveyoungsmodulus;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveshearmodulus;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveslidingfriction;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiverollingfriction;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveradius;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectivemass;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> modelparameterbeta;
=======

    void addParticleProperties(const std::shared_ptr<ParticlePropertyManager> manager);
    // Method to compute force between two spherical particles
    void computeSphereSphereForce(std::shared_ptr<SphereParticle> sphere1, std::shared_ptr<SphereParticle> sphere2,double timeStep);
    // Method to compute force between a sphere particle and a plane wall
    void computePlaneWallSphereForce(std::shared_ptr<PlaneWall> planeWall, std::shared_ptr<SphereParticle> sphere,double timeStep);
    
    void updateContactInformation();


private:

    std::unordered_map<int,std::unordered_map<int, ContactInformation>> contactinformationlist;
    std::unordered_map<int,std::unordered_map<int, ContactInformation>> nextcontactinformationlist;


    std::unordered_map<int,std::unordered_map<int, ContactInformation>> wallcontactinformationlist;
    std::unordered_map<int,std::unordered_map<int, ContactInformation>> nextwallcontactinformationlist;


    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveyoungsmodulus;

	std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveshearmodulus;

	std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveslidingfriction;

    std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiverollingfriction;

	std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectiveradius;

	std::map<PropertyTypeID, std::map<PropertyTypeID, double>> effectivemass;

	std::map<PropertyTypeID, std::map<PropertyTypeID, double>> modelparameterbeta;


>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
};