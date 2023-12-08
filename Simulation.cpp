// Simulation.cpp
#include "Simulation.h"

Simulation::Simulation(int numSpheres, double sphere_radius, double sphere_wall_radius, double Ec, double mass, double e)
    : mass(mass), numSpheres(numSpheres),
      contactforce(Ec, mass, e, sphere_radius), 
      spatialgrid(2 * sphere_wall_radius, 2 * sphere_wall_radius, 2 * sphere_wall_radius, 2.2 * sphere_radius)
{
    spherewall.radius = sphere_wall_radius;
    spherewall.x = sphere_wall_radius;
    spherewall.y = sphere_wall_radius;
    spherewall.z = sphere_wall_radius;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 2 * spherewall.radius);
    for (int i = 0; i < numSpheres; ++i)
    {
        Sphere sphere;
        bool validPosition = false;
        do
        {
            sphere.x = dis(gen);
            sphere.y = dis(gen);
            sphere.z = dis(gen);
            sphere.radius = sphere_radius;

            validPosition = isInside(sphere, spherewall);
            for (const auto &existsphere : spheres)
            {
                if (isOverlapping(sphere, existsphere))
                {
                    validPosition = false;
                    break;
                }
            }

        } while (!validPosition);
        sphere.id = i;
        sphere.vx = 0;
        sphere.vy = 0;
        sphere.vz = 0;

        sphere.wx = 0;
        sphere.wy = 0;
        sphere.wz = 0;

        sphere.forcex = 0;
        sphere.forcey = 0;
        sphere.forcez = 0;

        sphere.torquex = 0;
        sphere.torquey = 0;
        sphere.torquez = 0;

        spheres.push_back(sphere);
    }
    dt = 0.00001;
}

void Simulation::Update()
{
    spherewall.radius -= 0.1*dt;
    std::unordered_map<int, std::vector<int>> contact_paris;
    spatialgrid.broad_check(spheres, contact_paris);

    

    for (auto &contact_list : contact_paris)
    {
        int Main_sphere_id = contact_list.first;
        for (auto &next_sphere_id : contact_list.second)
        {
            double distance = sqrt((spheres[Main_sphere_id].x - spheres[next_sphere_id].x) * (spheres[Main_sphere_id].x - spheres[next_sphere_id].x) + (spheres[Main_sphere_id].y - spheres[next_sphere_id].y) * (spheres[Main_sphere_id].y - spheres[next_sphere_id].y) + (spheres[Main_sphere_id].z - spheres[next_sphere_id].z) * (spheres[Main_sphere_id].z - spheres[next_sphere_id].z));
            double overlap = (spheres[Main_sphere_id].radius + spheres[next_sphere_id].radius) - distance;

            if (overlap > 0)
            {
                contactforce.compute(spheres[Main_sphere_id], spheres[next_sphere_id], overlap);
            }
        }
    }

    for (auto &sphere : spheres)
    {
        double distance = sqrt((sphere.x - spherewall.x) * (sphere.x - spherewall.x) 
        + (sphere.y - spherewall.y) * (sphere.y - spherewall.y)
         + (sphere.z - spherewall.z) * (sphere.z - spherewall.z));
        double overlap = distance - (spherewall.radius - sphere.radius);
        if(overlap > 0)
        {
            contactforce.compute(sphere, spherewall, overlap);
        }
    }

    for (auto &sphere : spheres)
    {
        sphere.vx += sphere.forcex / mass * dt;
        sphere.vy += sphere.forcey / mass * dt;
        sphere.vz += sphere.forcez / mass * dt;

        sphere.x += sphere.vx * dt;
        sphere.y += sphere.vy * dt;
        sphere.z += sphere.vz * dt;

        sphere.forcex = 0;
        sphere.forcey = 0;
        sphere.forcez = 0;


    }
}

const std::vector<Sphere> &Simulation::GetSpheres() const
{
    return spheres;
}

const SphereWall &Simulation::GetSphereWall() const
{
    return spherewall;
}

bool Simulation::isInside(const Sphere &sphere, const SphereWall &spherewall)
{
    double distance = sqrt(pow(sphere.x - spherewall.x, 2) +
                           pow(sphere.y - spherewall.y, 2) +
                           pow(sphere.z - spherewall.z, 2));
    return distance + sphere.radius < spherewall.radius;
}

bool Simulation::isOverlapping(const Sphere &sphere1, const Sphere &sphere2)
{
    double distance = sqrt(pow(sphere2.x - sphere1.x, 2) +
                           pow(sphere2.y - sphere1.y, 2) +
                           pow(sphere2.z - sphere1.z, 2));
    return distance < (sphere1.radius + sphere2.radius);
}
