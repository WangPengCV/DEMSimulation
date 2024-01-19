#pragma once
#include "ParticlePropertyManager.h"
#include "PlaneWall.h"
#include "SphereParticle.h"
#include "SphereCylinderBond.h"
#include "Fiber.h"
#include "GridBasedContactDetection.h"
#include "ContactForce.h"
#include "RectangularContainer.h"
#include "SpringOscillator.h"
#include "FreeMotionTask.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <random>

class DEMProperties
{

public:
    DEMProperties();

    void loadFromFile(const std::string &filename);

    void saveToFile(const std::string &filename) const;

    std::shared_ptr<ParticlePropertyManager> getParticleManager() const
    {
        return manager;
    }

    std::shared_ptr<ContactForce> getContactForce() const
    {
        return contactforce;
    }

    const std::vector<std::shared_ptr<PlaneWall>> &getPlaneWall() const
    {
        return planewalls;
    }

    const std::vector<std::shared_ptr<SphereParticle>> &getsphereParticles() const
    {
        return sphereparticles;
    }

    const std::vector<std::shared_ptr<SphereParticle>> &getfibersphereParticles() const
    {
        return fibershpereparticles;
    }
    const std::vector<std::shared_ptr<SphereCylinderBond>> &getfiberbonds() const
    {
        return fiberbonds;
    }
    double getTotalTime() const
    {
        return totalTime;
    }

    double getTimestep() const
    {
        return timestep;
    }

    int getShowInterval() const
    {
        return showInterval;
    }

    double getAverageVelocity();

    const std::map<int, Eigen::Vector3d> &getSpecificForces() const
    {
        return specificforces;
    }

    const Eigen::Vector3d &getGlobalForces() const
    {
        return globalforces;
    }
    const Eigen::Vector3d &getSimulationDimensions() const
    {
        return simulationdimensions;
    }

    const std::shared_ptr<RectangularContainer> &getRectangularContainer() const
    {
        return rectangularcontainer;
    }
    bool isGenerateComplete() const { return (gernerateSphereFlag && gernerateFiberFlag) ; }
    // function for DEM
    void initialSimulation();
    void generateRemainingParticles();

    void handleCollisions();
    void applyExternalForces();
    void motion();

    void initial_task();
    // task function
    void dotask();
    

private:
    bool parseVector3d(std::istringstream &iss, Eigen::Vector3d &vec);

    void parseLine(const std::string &line);

    void line_process(std::string &line);

    void parseSphereProperties(std::istringstream &iss);

    void parsePlanewallProperties(std::istringstream &iss);

    void parseFiberProperties(std::istringstream &iss);

    void parsePlaneWall(std::istringstream &iss);

    void parseRectangularContainer(std::istringstream &iss);

    void parseRandomParticle(std::istringstream &iss);

    void parseSpecificParticle(std::istringstream &iss);

    std::shared_ptr<ParticlePropertyManager> manager;

    std::shared_ptr<ContactForce> contactforce;

    std::vector<std::shared_ptr<PlaneWall>> planewalls;

    std::shared_ptr<RectangularContainer> rectangularcontainer;

    std::vector<std::shared_ptr<SphereParticle>> sphereparticles;

    std::vector<std::shared_ptr<Fiber>> fibers; 

    std::vector<std::shared_ptr<SphereCylinderBond>> fiberbonds;

    std::vector<std::shared_ptr<SphereParticle>> fibershpereparticles;

    std::shared_ptr<GridBasedContactDetection> gridbasedcontactdetection;


    double timestep;

    double totalTime;

    int showInterval;

    Eigen::Vector3d gravity;

    std::map<int, Eigen::Vector3d> specificforces;

    Eigen::Vector3d globalforces;

    Eigen::Vector3d simulationdimensions;

    bool gernerateSphereFlag;

    bool gernerateFiberFlag;


    std::map<PropertyTypeID, std::string> gernerateInfor;

    std::unordered_map<std::string,std::unordered_map<int,std::shared_ptr<DynamicBoundaryMover>>> dynamicBoundaryTask;

    std::unordered_map<std::string,std::unordered_map<int,std::shared_ptr<EventTask>>> eventTasks;




};
