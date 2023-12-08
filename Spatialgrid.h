#pragma once
#include <unordered_map>
#include <unordered_set>
#include "Datatype.h"
#include <vector>
class Spatialgrid
{
public:

    Spatialgrid(double domain_x,double domain_y,double domain_z,double cellsize);

    void broad_check(const std::vector<Sphere>& spheres,std::unordered_map<int,std::vector<int>>& contact_paris);


private: 

    void cellmapping(const std::vector<Sphere>& spheres);

    int getCellIndex(double x, double y, double z);

    double domain_x;
    double domain_z;
    double domain_y;

    double cellsize_x;
    double cellsize_y;
    double cellsize_z;

    int cellnumber_x;
    int cellnumber_y;
    int cellnumber_z;

    std::unordered_map<int, std::unordered_set<int>> cells;


};

