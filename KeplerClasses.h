#ifndef KEPLERCLASSES_H
#define KEPLERCLASSES_H

#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace arma;


// ============================================== CLASS: CELESTIALBODY ============================================== //
class CelestialBody{
public:
    vec position;
    vec velocity;
    double mass;

    CelestialBody();
    CelestialBody(vec position, vec velocity, double mass);
};  // End of CelestialBody class declaration


// ============================================== CLASS: SYSTEM ===================================================== //
class System {
public:
    int     nBodies;
    double  dt;
    CelestialBody bodies[200];

    // Initialization
    System();
    void addBody(vec position, vec velocity, double mass);
    void definingBodies();
    vec  defineSystem;
    vec  energyAngMom;

    // Solving functions
    void solver(double nYears);
    void RK4();
    vec  diffEq(vec k);
    void calculateForcesAndEnergy(vec posVel, vec &Fx, vec &Fy, vec &Fz);

};  // End of System class declariation

#endif // KEPLERCLASSES_H
