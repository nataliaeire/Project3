#ifndef SYSTEM_H
#define SYSTEM_H
#include <fstream>
#include <vector>
#include <vec3.h>
#include <celestialbody.h>

class System
{
private:
    void    setG(bool cluster);
    void    actuallyCalculatingForces(CelestialBody &body, int n);

public:
    std::vector<CelestialBody> bodies;
    std::vector<CelestialBody*> bodies1;
    std::vector<CelestialBody*> bodies2;
    std::vector<CelestialBody*> bodies3;
    std::vector<CelestialBody*> bodies4;

    // Adding system
    System();
    void    addBody(vec3 position, vec3 velocity, double mass);
    void    addBody(double x, double y, double z, double vx, double vy, double vz, double mass);
    void    addSystem(std::fstream &file);
    void    addRandomSystem(int numberOfObjects, int sphereRadius);

    // Initialization of system
    void    conserveMomentum();
    void    sortBodiesIntoGroups();

    // Qualities of system
    int     numberOfBodies();
    double  G;
    double  totalMass;
    double  density;
    vec3    angularMomentum;
    vec3    momentum;
    double  totalEnergy();
    double  potentialEnergy;
    double  kineticEnergy;

    // Calculating qualities of the system
    //void    setG(bool cluster);
    void    calculateForcesAndEnergy();
    void    calculateForcesUsingGR();
    void    calculateForcesAdaptively(int n);  
};  // End of System class declariation

#endif // SYSTEM_H
