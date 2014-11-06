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

public:
    std::vector<CelestialBody> bodies;
    std::vector<CelestialBody*> bodies1;
    std::vector<CelestialBody*> bodies2;
    std::vector<CelestialBody*> bodies3;
    std::vector<CelestialBody*> bodies4;
    bool    smoothing;

    // Adding system
    System();
    void    addBody(vec3 position, vec3 velocity, double mass);
    void    addBody(double x, double y, double z, double vx, double vy, double vz, double mass);
    void    addSystem(std::fstream &file);
    void    addRandomSystem(int numberOfObjects, double sphereRadius);

    // Initialization of system
    void    conserveMomentum();
    void    resetEnergy();
    void    sortBodiesIntoGroups();

    // Qualities of system
    int     numberOfBodies();
    int     numberOfBoundBodies();
    double  G;
    double  totalMass;
    double  density;
    vec3    angularMomentum;
    vec3    momentum;
    double  totalEnergy();
    void    virialEnergy();
    double  potentialEnergy;
    double  kineticEnergy;
    double  boundKineticEnergy;
    double  boundPotentialEnergy;
    double  virialKineticEnergy;
    double  virialPotentialEnergy;

    // Calculating qualities of the system
/*    double  densityAsAFunctionOfRadius(double radius);
    double  deviationAndAverageDistanceBound();         // Not yet implemented */
    void    calculateForcesAndEnergy();
    void    calculateForcesUsingGR();
    void    actuallyCalculatingForces(CelestialBody &body, int n);
};  // End of System class declariation

#endif // SYSTEM_H
