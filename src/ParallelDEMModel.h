#include <mpi.h>
#include "DEMModel.h"
#include <filesystem>
#include <iostream>

class ParallelDEMModel {

public:struct ParticleData {
    int id;
    int state;
    double position[3];
    double velocity[3];
    double omega[3];
    double force[3];
    double torque[3];
};
    ParallelDEMModel(const std::string& filename);
    void runSimulation();
private:
    void domainDecomposition();
    std::shared_ptr<DEMProperties> DEMproperties;
    std::shared_ptr<Visualization> vis;
    int world_rank, world_size;
};