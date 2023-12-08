#include "Contactforce.h"

Contactforce::Contactforce(double Ec, double mass, double e, double r)
{
    effective_mass = mass / 2;
    effective_r = r / 2;
    effective_contact_modulus = Ec / (2 * (1 - 0.3 * 0.3));
    beta = log(e) / sqrt(3.14 * 3.14 + log(e) * log(e));
}

void Contactforce::compute(Sphere &sphere1, Sphere &sphere2, double overlap)
{
    Eigen::Vector3d normal = Eigen::Vector3d(sphere2.x - sphere1.x, sphere2.y - sphere1.y, sphere2.z - sphere1.z);

    normal.normalize();

    Eigen::Vector3d v1 = Eigen::Vector3d(sphere1.vx, sphere1.vy, sphere1.vz);

    Eigen::Vector3d v2 = Eigen::Vector3d(sphere2.vx, sphere2.vy, sphere2.vz);

    Eigen::Vector3d w1 = Eigen::Vector3d(sphere1.wx, sphere1.wy, sphere1.wz);

    Eigen::Vector3d w2 = Eigen::Vector3d(sphere2.wx, sphere2.wy, sphere2.wz);

    Eigen::Vector3d relative_v = (v1 - v2) + (sphere1.radius * w1 + sphere2.radius * w2).cross(normal);

    Eigen::Vector3d normal_v = relative_v.dot(normal) * normal;

    double kn = 1.3333 * effective_contact_modulus * sqrt(effective_r * overlap);

    double sn = 2 * effective_contact_modulus * sqrt(effective_r * overlap);

    double kdamping = -1.8257 * beta * sqrt(sn * effective_mass);

    Eigen::Vector3d normal_contact_force = kn * overlap * normal;

    Eigen::Vector3d normal_damping_force = kdamping * normal_v;

    Eigen::Vector3d normal_force = normal_contact_force + normal_damping_force;

    sphere1.forcex -= normal_force.x();
    sphere1.forcey -= normal_force.y();
    sphere1.forcez -= normal_force.z();

    sphere2.forcex += normal_force.x();
    sphere2.forcey += normal_force.y();
    sphere2.forcez += normal_force.z();

}


void Contactforce::compute(Sphere& sphere,SphereWall& spherewall,double overlap)
{
    Eigen::Vector3d normal = Eigen::Vector3d(spherewall.x - sphere.x, spherewall.y - sphere.y, spherewall.z - sphere.z);

    normal.normalize();

    Eigen::Vector3d v1 = Eigen::Vector3d(sphere.vx, sphere.vy, sphere.vz);

    Eigen::Vector3d w1 = Eigen::Vector3d(sphere.wx, sphere.wy, sphere.wz);

    Eigen::Vector3d relative_v = v1  - (sphere.radius * w1).cross(normal);

    Eigen::Vector3d normal_v = relative_v.dot(normal) * normal;

    double kn = 1.3333 * effective_contact_modulus * sqrt(effective_r * overlap);

    double sn = 2 * effective_contact_modulus * sqrt(effective_r * overlap);

    double kdamping = -1.8257 * beta * sqrt(sn * effective_mass);

    Eigen::Vector3d normal_contact_force = kn * overlap * normal;

    Eigen::Vector3d normal_damping_force = kdamping * normal_v;

    Eigen::Vector3d normal_force = normal_contact_force + normal_damping_force;

    sphere.forcex += normal_force.x();
    sphere.forcey += normal_force.y();
    sphere.forcez += normal_force.z();

    
}
