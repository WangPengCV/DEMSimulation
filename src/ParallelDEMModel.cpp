#include "ParallelDEMModel.h"
#include <mpi.h>
#include <filesystem>
#include <iostream>

ParallelDEMModel::ParallelDEMModel(const std::string &filename) {
    // 获取当前进程的等级和总进程数
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // 只有主进程负责读取文件和初始化
    if (world_rank == 0) {
        DEMproperties = std::make_shared<DEMProperties>();
        DEMproperties->loadFromFile(filename);
        vis = std::make_shared<Visualization>(*DEMproperties);
    }
    // 这里可能需要将DEMproperties的一些数据分发给其他进程
}

void ParallelDEMModel::runSimulation() {
    // 确保所有进程同步至此点
    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == 0) {
        DEMproperties->initialSimulation();
    }

    // 广播或分发必要的模拟参数到其他进程
    // 注意：你需要自己实现这部分逻辑，这可能包括DEMproperties中数据的序列化和广播

    double currentTime = 0.0;
    double timeStep; // 这将从主进程获取
    double totalTime; // 这将从主进程获取

    // 从主进程获取timeStep和totalTime
    // MPI_Bcast(&timeStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // MPI_Bcast(&totalTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    while (currentTime < totalTime) {
        // 这里插入并行处理逻辑和数据同步代码

        currentTime += timeStep;
    }

    // 主进程可能需要收集所有进程的结果
    if (world_rank == 0) {
        // 收集数据，进行最终的处理
    }
}
