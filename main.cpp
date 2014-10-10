#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <system.h>
#include <integrator.h>
#include <printing.h>

using namespace std;
using namespace arma;

int main()
{
    // Integration specifications
    double dt       = 0.01;
    double T        = 100;
    double nSteps   = T/dt;

    // Creating file to write to
    ofstream positionOut;
    ofstream velocityOut;
    ofstream energyAngMomOut;
    positionOut.open("position.dat");
    velocityOut.open("velocity.dat");
    energyAngMomOut.open("energyangmom.dat");

    // Qualities of the system we will be exploring
    System solarSystem;
    solarSystem.addBody(1, 0, 0, 0, 2*M_PI, 0, 3e-6);
    solarSystem.addBody(0, 0, 0, 0, -M_PI*6e-6, 0, 1);

    Integrator solver;
    CelestialBody &earth = solarSystem.bodies[0];

    cout << "Earth position before simulation: ["
         << earth.position[0] << ", "
         << earth.position[1] << ", "
         << earth.position[2] << "]" << endl;

    Printing printer("balle");
    printer.printingPosition(solarSystem);

    // Performing RK4 on the system
    for(int i = 0; i < nSteps; i++){
        solver.RK4(solarSystem, dt);
        printer.printingPosition(solarSystem);
    }

    cout << "Earth position after simulation: ["
         << earth.position[0] << ", "
         << earth.position[1] << ", "
         << earth.position[2] << "]" << endl;

    return 0;
}
