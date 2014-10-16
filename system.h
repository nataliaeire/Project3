#ifndef SYSTEM_H
#define SYSTEM_H
#include <fstream>
#include <vector>
#include <vec3.h>
#include <celestialbody.h>

class System
{
public:
    std::vector<CelestialBody> bodies;

    // Initialization
    System();
    void    addBody(vec3 position, vec3 velocity, double mass);
    void    addBody(double x, double y, double z, double vx, double vy, double vz, double mass);
    void    addSystem(std::fstream &file);

    // Calculated qualities of the system
    void    conserveMomentum();
    void    calculateForcesAndEnergy();
    void    calculateForcesUsingGR();
    int     numberOfBodies();
    double  totalEnergy();
    double  potentialEnergy;
    double  kineticEnergy;
    vec3    angularMomentum;
    vec3    momentum;
};  // End of System class declariation

#endif // SYSTEM_H
