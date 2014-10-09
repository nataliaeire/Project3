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
    int nBodies;
    CelestialBody bodies[200];

    // Initialization
    System();
    void addBody(vec position, vec velocity, double mass);
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
