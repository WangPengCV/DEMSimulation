#pragma once
#include "cmath"
#include "Datatype.h"
#include <Eigen/Dense>

class Contactforce
{
public:

    Contactforce(double Ec, double mass, double e, double r);

    void compute(Sphere& sphere1,Sphere& sphere2,double overlap);

    void compute(Sphere& sphere,SphereWall& spherewall,double overlap);


private:

    double effective_mass;

    double effective_contact_modulus;

    double effective_r;

    double  beta;

    

};