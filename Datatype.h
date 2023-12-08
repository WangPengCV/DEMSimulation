#pragma once
struct Sphere
{
    int id;
    double x;
    double y;
    double z;

    double vx;
    double vy;
    double vz;

    double wx;
    double wy;
    double wz;

    double forcex;
    double forcey;
    double forcez;

    double torquex;
    double torquey;
    double torquez;

    double radius;
};
struct SphereWall
{
    double x;
    double y;
    double z;
    double radius;
};