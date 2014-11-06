#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H
#include "vec3.h"



class CelestialBody
{
public:
    int     index;
    vec3    position;
    vec3    velocity;
    vec3    force;
    double  mass;
    double  PE;
    double  KE;
    vec3    angMom;

    // Initialisation
    CelestialBody(int index, vec3 position, vec3 velocity, double mass);
    void    resetForce();
    void    resetEnergy();
    vec3    acceleration();
    double  totalEnergyOfBody();

};  // End of CelestialBody class declaration

#endif // CELESTIALBODY_H
