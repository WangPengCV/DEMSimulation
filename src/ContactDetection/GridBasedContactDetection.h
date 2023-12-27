#pragma once
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <SphereParticle.h>
class GridBasedContactDetection
{
public:

    GridBasedContactDetection();

    void initial(double domain_x,double domain_y,double domain_z,double gridsize);

    void ParticleBroadPhase(const std::vector<std::shared_ptr<Particle>> &particles,std::unordered_map<int,std::unordered_set<int>>& contact_paris);

    void planewallBroadPhase(const std::vector<std::shared_ptr<Particle>> &particles, std::vector<std::shared_ptr<PlaneWall>>& planewalls,
                             std::unordered_map<int,std::unordered_set<int>>& contact_paris);

    double getDomainX  () const
    {
        return domainX;
    }

    double getDomainY  () const
    {
        return domainY;
    }

    double getDomainZ  () const
    {
        return domainZ;
    }

    double getGridSizeX() const
    {
        return gridSizeX;
    }

    double getGridSizeY() const
    {
        return gridSizeY;
    }

    double getGridSizeZ() const
    {
        return gridSizeZ;
    }

    int getNumberOfGridX() const
    {
        return numberOfGridX;
    }

    int getNumberOfGridY() const
    {
        return numberOfGridY;
    }

    int getNumberOfGridZ() const
    {
        return numberOfGridZ;
    }

private: 

    void assignParticleToGrid(const std::vector<std::shared_ptr<Particle>> &particles);

    int getGridIndex(double x, double y, double z);

    double domainX;
    double domainY;
    double domainZ;

    double gridSizeX;
    double gridSizeY;
    double gridSizeZ;
    double gridsize;
    int numberOfGridX;
    int numberOfGridY;
    int numberOfGridZ;

    std::unordered_map<int, std::unordered_set<int>> gridSaveParticles;

    std::unordered_map<int, std::unordered_set<int>> planeSaveGrid;



};

