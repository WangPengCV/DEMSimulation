#include "Spatialgrid.h"

Spatialgrid::Spatialgrid(double domain_x, double domain_y, double domain_z, double cellsize)
    : domain_x(domain_x), domain_y(domain_y), domain_z(domain_z)
{
    cellnumber_x = static_cast<int>(std::ceil(domain_x / cellsize));
    cellnumber_y = static_cast<int>(std::ceil(domain_y / cellsize));
    cellnumber_z = static_cast<int>(std::ceil(domain_z / cellsize));

    cellsize_x = domain_x / cellnumber_x;
    cellsize_y = domain_y / cellnumber_y;
    cellsize_z = domain_z / cellnumber_z;
}

int Spatialgrid::getCellIndex(double x, double y, double z)
{
    int cellX = static_cast<int>(x / cellsize_x);
    int cellY = static_cast<int>(y / cellsize_y);
    int cellZ = static_cast<int>(z / cellsize_z);

    return cellX + cellZ * cellnumber_x + cellY * cellnumber_x * cellnumber_z;
}

void Spatialgrid::cellmapping(const std::vector<Sphere> &spheres)
{
    for (const auto &sphere : spheres)
    {
        cells[getCellIndex(sphere.x - sphere.radius, sphere.y - sphere.radius, sphere.z - sphere.radius)].insert(sphere.id);
        cells[getCellIndex(sphere.x - sphere.radius, sphere.y - sphere.radius, sphere.z + sphere.radius)].insert(sphere.id);
        cells[getCellIndex(sphere.x + sphere.radius, sphere.y - sphere.radius, sphere.z - sphere.radius)].insert(sphere.id);
        cells[getCellIndex(sphere.x + sphere.radius, sphere.y - sphere.radius, sphere.z + sphere.radius)].insert(sphere.id);

        cells[getCellIndex(sphere.x - sphere.radius, sphere.y + sphere.radius, sphere.z - sphere.radius)].insert(sphere.id);
        cells[getCellIndex(sphere.x - sphere.radius, sphere.y + sphere.radius, sphere.z + sphere.radius)].insert(sphere.id);
        cells[getCellIndex(sphere.x + sphere.radius, sphere.y + sphere.radius, sphere.z - sphere.radius)].insert(sphere.id);
        cells[getCellIndex(sphere.x + sphere.radius, sphere.y + sphere.radius, sphere.z + sphere.radius)].insert(sphere.id);
    }
}

void Spatialgrid::broad_check(const std::vector<Sphere> &spheres,std::unordered_map<int,std::vector<int>>& contact_paris)
{
    cellmapping(spheres);

    for (auto &cell : cells)
    {
        const auto &sphere_ids = cell.second;
        for (const auto& main_id : sphere_ids)
        {
            for(const auto next_id : sphere_ids)
            {
                if(main_id < next_id)
                {
                    contact_paris[main_id].push_back(next_id);
                }

            }
        }
    }
}
