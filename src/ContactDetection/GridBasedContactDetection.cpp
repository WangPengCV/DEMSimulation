#include "GridBasedContactDetection.h"

GridBasedContactDetection::GridBasedContactDetection()
{
}
void GridBasedContactDetection::initial(double domain_x, double domain_y, double domain_z, double girdsize)
{
    domainX = domain_x;
    domainY = domain_y;
    domainZ = domain_z;
    numberOfGridX = static_cast<int>(std::ceil(domainX / girdsize));
    numberOfGridY = static_cast<int>(std::ceil(domainY / girdsize));
    numberOfGridZ = static_cast<int>(std::ceil(domainZ / girdsize));

    gridSizeX = domainX / numberOfGridX;
    gridSizeY = domainX / numberOfGridY;
    gridSizeZ = domainZ / numberOfGridZ;

    gridsize = std::min(std::min(gridSizeX, gridSizeY), gridSizeZ);
}

int GridBasedContactDetection::getGridIndex(double x, double y, double z)
{
    int cellX = static_cast<int>(x / gridSizeX);
    int cellY = static_cast<int>(y / gridSizeY);
    int cellZ = static_cast<int>(z / gridSizeZ);

    return cellX + cellZ * numberOfGridX + cellY * numberOfGridX * numberOfGridZ;
}

void GridBasedContactDetection::assignParticleToGrid(const std::vector<std::shared_ptr<SphereParticle>> &sphereparticles,
                                                     const std::vector<std::shared_ptr<SphereCylinderBond>> &fiberbonds,
                                                     const std::vector<std::shared_ptr<SphereParticle>> &fibersphereparticles)
{
    gridSaveSphereParticles.clear(); // Clear previous data
    gridSaveFiberBonds.clear();
    for (const auto &particle : sphereparticles)
    {
        int category = particle->getType().getCategory();
        auto manager = particle->getParticlePropertyManager();

        auto radius = manager->getSphereProperties(particle->getType())->getRadius();
        int minCellX = static_cast<int>((particle->getPosition().x() - radius) / gridSizeX);
        int minCellY = static_cast<int>((particle->getPosition().y() - radius) / gridSizeY);
        int minCellZ = static_cast<int>((particle->getPosition().z() - radius) / gridSizeZ);

        int maxCellX = static_cast<int>((particle->getPosition().x() + radius) / gridSizeX);
        int maxCellY = static_cast<int>((particle->getPosition().y() + radius) / gridSizeY);
        int maxCellZ = static_cast<int>((particle->getPosition().z() + radius) / gridSizeZ);

        // Iterate over the range of grid cells and assign the particle's ID
        for (int x = minCellX; x <= maxCellX; x++)
        {
            for (int z = minCellZ; z <= maxCellZ; z++)
            {
                for (int y = minCellY; y <= maxCellY; y++)
                {

                    int gridIndex = x + z * numberOfGridX + y * numberOfGridX * numberOfGridZ;
                    gridSaveSphereParticles[gridIndex].insert(particle->getId());
                }
            }
        }
    }
    for (const auto &particle : fiberbonds)
    {
        int category = particle->getType().getCategory();
        auto manager = particle->getParticlePropertyManager();

        auto radius = manager->getFiberProperties(particle->getType())->getRadius();

        int node1Id = particle->getNode1();
        int node2Id = particle->getNode2();

        Eigen::Vector3d startposition = fibersphereparticles[node1Id]->getPosition();
        Eigen::Vector3d endposition = fibersphereparticles[node2Id]->getPosition();

        int numSegments = static_cast<int>(std::ceil(manager->getFiberProperties(particle->getType())->getElementlength() / (2 * radius)));
        Eigen::Vector3d segmentVector = (endposition - startposition) / numSegments;
        for (int i = 0; i <= numSegments; ++i)
        {
            Eigen::Vector3d center;
            if (i < numSegments)
            {
                center = startposition + i * segmentVector;
            }
            else
            {
                // For the last segment, set the center to the endposition
                center = endposition;
            }

            int minCellX = static_cast<int>((center.x() - radius) / gridSizeX);
            int minCellY = static_cast<int>((center.y() - radius) / gridSizeY);
            int minCellZ = static_cast<int>((center.z() - radius) / gridSizeZ);

            int maxCellX = static_cast<int>((center.x() + radius) / gridSizeX);
            int maxCellY = static_cast<int>((center.y() + radius) / gridSizeY);
            int maxCellZ = static_cast<int>((center.z() + radius) / gridSizeZ);

            // Iterate over the range of grid cells and assign the particle's ID
            for (int x = minCellX; x <= maxCellX; x++)
            {
                for (int z = minCellZ; z <= maxCellZ; z++)
                {
                    for (int y = minCellY; y <= maxCellY; y++)
                    {

                        int gridIndex = x + z * numberOfGridX + y * numberOfGridX * numberOfGridZ;
                        gridSaveFiberBonds[gridIndex].insert(particle->getId());
                    }
                }
            }
        }
    }
}

void GridBasedContactDetection::ParticleBroadPhase(const std::vector<std::shared_ptr<SphereParticle>> &sphereparticles,
                                                   const std::vector<std::shared_ptr<SphereCylinderBond>> &fiberbonds,
                                                   const std::vector<std::shared_ptr<SphereParticle>> &fibersphereparticles,
                                                   std::unordered_map<int, std::unordered_set<int>> &SScontactparis,
                                                   std::unordered_map<int, std::unordered_set<int>> &SFcontactparis,
                                                   std::unordered_map<int, std::unordered_set<int>> &FFcontactparis)
{
    assignParticleToGrid(sphereparticles, fiberbonds, fibersphereparticles);
    for (const auto &grid_entry : gridSaveSphereParticles)
    {
        const auto &sphere_in_grid = grid_entry.second;

        for (auto it1 = sphere_in_grid.begin(); it1 != sphere_in_grid.end(); ++it1)
        {

            for (auto it2 = std::next(it1); it2 != sphere_in_grid.end(); ++it2)
            {
                int id1 = *it1;
                int id2 = *it2;
                SScontactparis[std::min(id1, id2)].insert(std::max(id1, id2));
            }
            if (gridSaveFiberBonds.find(grid_entry.first) != gridSaveFiberBonds.end())
            {
                for (auto it3 = gridSaveFiberBonds[grid_entry.first].begin(); it3 != gridSaveFiberBonds[grid_entry.first].end(); ++it3)
                {
                    int sphere_id = *it1;
                    int fiber_id = *it3;
                    SFcontactparis[sphere_id].insert(fiber_id);
                }
            }
        }
    }
    for (const auto &grid_entry : gridSaveFiberBonds)
    {
        const auto &fibers_in_grid = grid_entry.second;

        for (auto it1 = fibers_in_grid.begin(); it1 != fibers_in_grid.end(); ++it1)
        {
            auto it2 = it1;
            ++it2; // Start the second iterator ahead of the first
            for (; it2 != fibers_in_grid.end(); ++it2)
            {
                int id1 = *it1;
                int id2 = *it2;
                if(fiberbonds[id1]->getFiberId() != fiberbonds[id2]->getFiberId())
                {
                    FFcontactparis[std::min(id1, id2)].insert(std::max(id1, id2));
                }
            }
        }
    }
}

void GridBasedContactDetection::planewallBroadPhase(std::vector<std::shared_ptr<PlaneWall>> &planewalls,
                                                    std::unordered_map<int, std::unordered_set<int>> &PScontactpairs,
                                                    std::unordered_map<int, std::unordered_set<int>> &PFcontactpairs)
{

    for (auto &plane : planewalls)
    {
        if (plane->getState() == 1 || planeSaveGrid.find(plane->getId()) == planeSaveGrid.end())
        {
            plane->generateMesh(gridsize);
            std::vector<Eigen::Vector3d> meshvertices = plane->getMeshVertices();
            planeSaveGrid[plane->getId()].clear();
            std::unordered_set<int> &gridIndices = planeSaveGrid[plane->getId()];

            for (auto &vertice : meshvertices)
            {

                gridIndices.insert(getGridIndex(vertice.x(), vertice.y(), vertice.z()));
            }
        }
    }

    for (const auto &grid_entry : planeSaveGrid)
    {
        auto wall_id = grid_entry.first;
        const auto &grid_indices = grid_entry.second;
        auto &pscontactSet = PScontactpairs[wall_id];
        auto &pfcontactSet = PFcontactpairs[wall_id];

        for (int grid_id : grid_indices)
        {
            if (gridSaveSphereParticles.find(grid_id) != gridSaveSphereParticles.end())
            {
                const auto &particlesInGrid = gridSaveSphereParticles[grid_id];

                pscontactSet.insert(particlesInGrid.begin(), particlesInGrid.end());
            }

            if (gridSaveFiberBonds.find(grid_id) != gridSaveFiberBonds.end())
            {
                const auto &particlesInGrid = gridSaveFiberBonds[grid_id];

                pfcontactSet.insert(particlesInGrid.begin(), particlesInGrid.end());
            }
        }
    }
}
void GridBasedContactDetection::RectangularContainerBroadPhase(std::shared_ptr<RectangularContainer> &rectangularcontainer,
                                                               std::unordered_map<int, std::unordered_set<int>> &PScontactpairs,
                                                               std::unordered_map<int, std::unordered_set<int>> &PFcontactpairs)
{

    for (auto &plane : rectangularcontainer->getPlaneWall())
    {
        // if the plane can move freely or rectangularcontainerSaveGrid have no gird
        if (plane->getState() == 1 || rectangularcontainerSaveGrid.find(plane->getId()) == rectangularcontainerSaveGrid.end())
        {
            plane->generateMesh(gridsize);
            std::vector<Eigen::Vector3d> meshvertices = plane->getMeshVertices();
            rectangularcontainerSaveGrid[plane->getId()].clear();
            std::unordered_set<int> &gridIndices = rectangularcontainerSaveGrid[plane->getId()];

            for (auto &vertice : meshvertices)
            {
                gridIndices.insert(getGridIndex(vertice.x(), vertice.y(), vertice.z()));
            }
        }
    }

    for (const auto &grid_entry : rectangularcontainerSaveGrid)
    {
        auto wall_id = grid_entry.first;
        const auto &grid_indices = grid_entry.second;
        auto &pscontactSet = PScontactpairs[wall_id];
        auto &pfcontactSet = PFcontactpairs[wall_id];


        for (int grid_id : grid_indices)
        {
            if (gridSaveSphereParticles.find(grid_id) != gridSaveSphereParticles.end())
            {
                const auto &particlesInGrid = gridSaveSphereParticles[grid_id];

                pscontactSet.insert(particlesInGrid.begin(), particlesInGrid.end());
            }

            if (gridSaveFiberBonds.find(grid_id) != gridSaveFiberBonds.end())
            {
                const auto &particlesInGrid = gridSaveFiberBonds[grid_id];

                pfcontactSet.insert(particlesInGrid.begin(), particlesInGrid.end());
            }
        }
    }
}