#pragma once
#include "ParticlePropertyManager.h"
#include "PlaneWall.h"
#include "SphereParticle.h"
#include "GridBasedContactDetection.h"
#include "ContactForce.h"
#include "RectangularContainer.h"
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

    const std::vector<std::shared_ptr<Particle>> &getParticles() const
    {
        return particles;
    }

    double getTotalTime() const
    {
        return totalTime;
    }

    double getTimestep() const
    {
        return timestep;
    }

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

    const std::shared_ptr<RectangularContainer>& getRectangularContainer() const
    {
        return rectangularcontainer;
    }

    //function for DEM
    void initialContactDetection();
    void handleCollisions();
    void applyExternalForces();
    void motion();




private:
    bool parseVector3d(std::istringstream &iss, Eigen::Vector3d &vec);

    void parseLine(const std::string &line);

    void line_process(std::string &line);

    void parseSphereProperties(std::istringstream &iss);

    void parsePlanewallProperties(std::istringstream &iss);

    void parsePlaneWall(std::istringstream &iss);

    void parseRectangularContainer(std::istringstream &iss);


    void parseRandomParticle(std::istringstream &iss);

    void parseSpecificParticle(std::istringstream &iss);

    void generateRemainingParticles();


    std::shared_ptr<ParticlePropertyManager> manager;

    std::shared_ptr<ContactForce> contactforce;

    std::vector<std::shared_ptr<PlaneWall>> planewalls;

    std::shared_ptr<RectangularContainer> rectangularcontainer;

    std::vector<std::shared_ptr<Particle>> particles;

    std::shared_ptr<GridBasedContactDetection> gridbasedcontactdetection;


    double timestep;

    double totalTime;

    Eigen::Vector3d gravity;

    std::map<int, Eigen::Vector3d> specificforces;

    Eigen::Vector3d globalforces;

    Eigen::Vector3d simulationdimensions;


    bool gernerateFlag;
    std::map<int,std::string> gernerateInfor;
};
