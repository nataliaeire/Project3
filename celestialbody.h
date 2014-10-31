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

    // Initialisation
    CelestialBody(int index, vec3 position, vec3 velocity, double mass);
    void resetForce();
    vec3 acceleration();

};  // End of CelestialBody class declaration

#endif // CELESTIALBODY_H
