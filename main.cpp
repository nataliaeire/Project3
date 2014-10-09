#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <system.h>
#include <rk4.h>

using namespace std;
using namespace arma;

int main() {
    double dt = 0.01;
    double T = 10;
    double nSteps = T / dt;
    System solarSystem;
    solarSystem.addBody(1, 0, 0, 0, 2*M_PI, 0, 3e-6);
    solarSystem.addBody(0, 0, 0, 0, -M_PI*6e-6, 0, 1);

    RK4 solver;
    CelestialBody &earth = solarSystem.bodies[0];
    cout << "Earth position before simulation: [" << earth.position[0] << ", " << earth.position[1] << ", " << earth.position[2] << "]" << endl;
    for(int i = 0; i < nSteps; i++) {
        solver.integrate(solarSystem, dt);
        // cout << "Total energy: " << solarSystem.totalEnergy() << endl;
    }
    cout << "Earth position after simulation: [" << earth.position[0] << ", " << earth.position[1] << ", " << earth.position[2] << "]" << endl;




    //CelestialBody earth(1, 0, 0, 2*M_PI, 0, 0, 3e-6);
//    System solarsyst;
//    vec pos(2);
//    pos[0] = 1;
//    pos[1] = 0;
//    vec vel(2);
//    vel[0] = 0;
//    vel[1] = 2*M_PI;
//    double mass = 3e-6;

//    solarsyst.addBody(pos,vel,mass);
//    solarsyst.definingBodies();

//    solarsyst.solver(1);

    return 0;
}
