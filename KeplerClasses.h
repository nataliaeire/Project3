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
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double mass;

    CelestialBody();
    CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass);
};  // End of CelestialBody class declaration


// ============================================== CLASS: SYSTEM ===================================================== //
class System {
public:
    int nBodies;
    CelestialBody* bodies;

    // Initialization
    System(int n);
    vec defineSystem;

    // Initialization functions
    void twoBodies();
    void threeBodies();
    void definingBodies();

    // Solving functions
    void RK4(double h, ofstream* posvel, ofstream* enmom);
    vec  diffEq(vec, ofstream*);

};  // End of System class declariation

#endif // KEPLERCLASSES_H
