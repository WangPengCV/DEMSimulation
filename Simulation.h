#pragma once
#include <vector>
#include <cstdlib>
#include <ctime>
#include <random>
#include <iostream>
#include "Spatialgrid.h"
#include "Contactforce.h"
#include "Datatype.h"


class Simulation
{
public:
    explicit Simulation(int numSpheres, double sphere_radius, double sphere_wall_radius, double Ec, double mass, double e);
    void Update();
    const std::vector<Sphere> &GetSpheres() const;
    const SphereWall &GetSphereWall() const;

private:
    int numSpheres;

    std::vector<Sphere> spheres;

    SphereWall spherewall;

    bool isInside(const Sphere &sphere, const SphereWall &spherewall);

    bool isOverlapping(const Sphere &sphere1, const Sphere &sphere2);

    Contactforce contactforce;

    Spatialgrid spatialgrid;
    
    double dt;

    double mass;
};
