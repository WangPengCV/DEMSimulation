#pragma once
#include "ParticlePropertyManager.h"
#include "PlaneWall.h"
#include "SphereParticle.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <random>


class DEMProperties 
{

public:

    DEMProperties();
    
    void loadFromFile(const std::string& filename);

    void saveToFile(const std::string& filename) const;


    std::shared_ptr<ParticlePropertyManager> getParticleManager() const 
    {
        return manger;
    }

    const std::vector<std::shared_ptr<PlaneWall>>& getPlaneWall() const{
        return planewalls;
    }

    const std::vector<std::shared_ptr<Particle>>& getParticles() const{
        return particles;
    }

    double getTotalTime() const {
        return totalTime;
    }

    double getTimestep() const {
        return timestep;
    }

private:

    bool parseVector3d(std::istringstream& iss, Eigen::Vector3d& vec);

    void parseLine(const std::string& line); 

    void line_process(std::string &line);

	void parseSphereProperties(std::istringstream &iss);

	void parsePlanewallProperties(std::istringstream &iss);

	void parsePlaneWall(std::istringstream &iss);

    void parseRandomParticle(std::istringstream &iss);

    void parseSpecificParticle(std::istringstream &iss);




    std::shared_ptr<ParticlePropertyManager> manger;

    std::vector<std::shared_ptr<PlaneWall>> planewalls;

    std::vector<std::shared_ptr<Particle>> particles;

    double timestep;

    double totalTime;
    
};
