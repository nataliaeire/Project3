#ifndef SYSTEM_H
#define SYSTEM_H
#include <vector>
#include <celestialbody.h>
#include <armadillo>

class System
{
public:
    std::vector<CelestialBody> bodies;

    // Initialization
    System();
    void        addBody(arma::vec position, arma::vec velocity, double mass);
    void        addBody(double x, double y, double z, double vx, double vy, double vz, double mass);

    // Calculated qualities of the system
    void        calculateForcesAndEnergy();
    int         numberOfBodies();
    double      totalEnergy();
    double      potentialEnergy;
    double      kineticEnergy;
    arma::vec   angularMomentum;

};  // End of System class declariation

#endif // SYSTEM_H
