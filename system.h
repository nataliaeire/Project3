#ifndef SYSTEM_H
#define SYSTEM_H
#include <vector>
#include <celestialbody.h>
#include <vec3.h>

class System
{
public:
    std::vector<CelestialBody> bodies;

    // Initialization
    System();
    void        addBody(vec3 position, vec3 velocity, double mass);
    void        addBody(double x, double y, double z, double vx, double vy, double vz, double mass);

    // Calculated qualities of the system
    void        calculateForcesAndEnergy();
    int         numberOfBodies();
    double      totalEnergy();
    double      potentialEnergy;
    double      kineticEnergy;
    vec3        angularMomentum;

};  // End of System class declariation

#endif // SYSTEM_H
