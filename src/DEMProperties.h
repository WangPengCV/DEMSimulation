#pragma once
#include "ParticlePropertyManager.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

const double PI = 3.1415926535897932384626;

class DEMProperties 
{

public:
    DEMProperties();
    
    void loadFromFile(const std::string& filename);

    std::shared_ptr<ParticlePropertyManager> getParticleManage() const 
    {
        return manger;
    }

private:

    void parseLine(const std::string& line); 

    std::shared_ptr<ParticlePropertyManager> manger;
};
