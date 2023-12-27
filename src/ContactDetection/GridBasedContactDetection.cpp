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

void GridBasedContactDetection::assignParticleToGrid(const std::vector<std::shared_ptr<Particle>> &particles)
{
    gridSaveParticles.clear(); // Clear previous data
    for (const auto &particle : particles)
    {
        int category = particle->getType().getCategory();
        auto manager = particle->getParticlePropertyManager();
        auto typemapping = manager->gettypeMapping();
        if (typemapping[category] == ParticleType::SPHERE)
        {
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
                        gridSaveParticles[gridIndex].insert(particle->getId());
                    }
                }
            }
        }
    }
}

void GridBasedContactDetection::ParticleBroadPhase(const std::vector<std::shared_ptr<Particle>> &particles, std::unordered_map<int, std::unordered_set<int>> &pp_contact_paris)
{
    assignParticleToGrid(particles);
    for (const auto &grid_entry : gridSaveParticles)
    {
        const auto &particles_in_grid = grid_entry.second;

        for (auto it1 = particles_in_grid.begin(); it1 != particles_in_grid.end(); ++it1)
        {

            for (auto it2 = std::next(it1); it2 != particles_in_grid.end(); ++it2)
            {
                int id1 = *it1;
                int id2 = *it2;
                pp_contact_paris[std::min(id1, id2)].insert(std::max(id1, id2));
            }
        }
    }
}

void GridBasedContactDetection::planewallBroadPhase(const std::vector<std::shared_ptr<Particle>> &particles, std::vector<std::shared_ptr<PlaneWall>> &planewalls,
                                                    std::unordered_map<int, std::unordered_set<int>> &pw_contact_paris)
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
        auto &contactSet = pw_contact_paris[wall_id];

        for (int grid_id : grid_indices)
        {
            if (gridSaveParticles.find(grid_id) != gridSaveParticles.end())
            {
                const auto &particlesInGrid = gridSaveParticles[grid_id];
                
                contactSet.insert(particlesInGrid.begin(), particlesInGrid.end());
                
               
               
            }
        }
    }
}