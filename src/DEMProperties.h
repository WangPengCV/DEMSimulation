#pragma once
#include "ParticlePropertyManager.h"
#include "PlaneWall.h"
#include "SphereParticle.h"
<<<<<<< HEAD
#include "SphereCylinderBond.h"
#include "Fiber.h"
#include "GridBasedContactDetection.h"
#include "ContactForce.h"
#include "RectangularContainer.h"
#include "SpringOscillator.h"
#include "FreeMotionTask.h"
=======
#include "GridBasedContactDetection.h"
#include "ContactForce.h"
#include "RectangularContainer.h"
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
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

<<<<<<< HEAD
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
=======
    const std::vector<std::shared_ptr<Particle>> &getParticles() const
    {
        return particles;
    }

>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    double getTotalTime() const
    {
        return totalTime;
    }

    double getTimestep() const
    {
        return timestep;
    }

<<<<<<< HEAD
    int getShowInterval() const
    {
        return showInterval;
    }

    double getAverageVelocity();

=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
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

<<<<<<< HEAD
    const std::shared_ptr<RectangularContainer> &getRectangularContainer() const
    {
        return rectangularcontainer;
    }
    bool isGenerateComplete() const { return (gernerateSphereFlag && gernerateFiberFlag) ; }
    // function for DEM
    void initialSimulation();
    void generateRemainingParticles();

=======
    const std::shared_ptr<RectangularContainer>& getRectangularContainer() const
    {
        return rectangularcontainer;
    }

    //function for DEM
    void initialContactDetection();
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    void handleCollisions();
    void applyExternalForces();
    void motion();

<<<<<<< HEAD
    void initial_task();
    // task function
    void dotask();
    
=======


>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200

private:
    bool parseVector3d(std::istringstream &iss, Eigen::Vector3d &vec);

    void parseLine(const std::string &line);

    void line_process(std::string &line);

    void parseSphereProperties(std::istringstream &iss);

    void parsePlanewallProperties(std::istringstream &iss);

<<<<<<< HEAD
    void parseFiberProperties(std::istringstream &iss);

=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    void parsePlaneWall(std::istringstream &iss);

    void parseRectangularContainer(std::istringstream &iss);

<<<<<<< HEAD
=======

>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    void parseRandomParticle(std::istringstream &iss);

    void parseSpecificParticle(std::istringstream &iss);

<<<<<<< HEAD
=======
    void generateRemainingParticles();


>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    std::shared_ptr<ParticlePropertyManager> manager;

    std::shared_ptr<ContactForce> contactforce;

    std::vector<std::shared_ptr<PlaneWall>> planewalls;

    std::shared_ptr<RectangularContainer> rectangularcontainer;

<<<<<<< HEAD
    std::vector<std::shared_ptr<SphereParticle>> sphereparticles;

    std::vector<std::shared_ptr<Fiber>> fibers; 

    std::vector<std::shared_ptr<SphereCylinderBond>> fiberbonds;

    std::vector<std::shared_ptr<SphereParticle>> fibershpereparticles;
=======
    std::vector<std::shared_ptr<Particle>> particles;
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200

    std::shared_ptr<GridBasedContactDetection> gridbasedcontactdetection;


    double timestep;

    double totalTime;

<<<<<<< HEAD
    int showInterval;

=======
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
    Eigen::Vector3d gravity;

    std::map<int, Eigen::Vector3d> specificforces;

    Eigen::Vector3d globalforces;

    Eigen::Vector3d simulationdimensions;

<<<<<<< HEAD
    bool gernerateSphereFlag;

    bool gernerateFiberFlag;


    std::map<PropertyTypeID, std::string> gernerateInfor;

    std::unordered_map<std::string,std::unordered_map<int,std::shared_ptr<DynamicBoundaryMover>>> dynamicBoundaryTask;

    std::unordered_map<std::string,std::unordered_map<int,std::shared_ptr<EventTask>>> eventTasks;




=======

    bool gernerateFlag;
    std::map<int,std::string> gernerateInfor;
>>>>>>> 686cbfa3ebadc1d4aba7bce443978911f7964200
};
