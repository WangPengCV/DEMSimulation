#pragma once
#include "ParticlePropertyManager.h"
#include "PlaneWall.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>


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

private:

    bool parseVector3d(std::istringstream& iss, Eigen::Vector3d& vec);

    void parseLine(const std::string& line); 

    void line_process(std::string &line);

	void parseSphereProperties(std::istringstream &iss);

	void parsePlanewallProperties(std::istringstream &iss);

	void parsePlaneWall(std::istringstream &iss);

    std::shared_ptr<ParticlePropertyManager> manger;

    std::vector<PlaneWall> planewalls;

    
};
