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
    double dt       = 0.001;
    double T        = 1;
    double nSteps   = T/dt;

    // Qualities of the system we will be exploring are read from file
    // Note that file directory has to be changed accordingly for every computer
    fstream file("/uio/hume/student-u81/natales/Project3/Project3/hei.txt",ios_base::in);

    // Intialisation
    System SunEarthSystem;                      // Preparing system
    Integrator solving;                         // Preparing for allowing the system to develop
    Printing printer("SunEarth");               // Preparing for printing details about system to file

    SunEarthSystem.addSystem(file);             // Creating system
    printer.printingAll(SunEarthSystem);        // Printing intitial values to file

    int counter = 1;                            // Counting parameter to print message to screen inside for-loop

    // Performing RK4 on the system
    for(int i = 0; i < nSteps; i++){
        solving.RK4(SunEarthSystem, dt);        // Solving the problem using the RK4-method
        printer.printingAll(SunEarthSystem);    // Printing everything to file
        counter ++;

        // Printing a message to screen to let the user know how far the program has come
        if(counter == floor(0.25*nSteps))   cout << "25 % of the integration is performed" << endl;
        if(counter == floor(0.50*nSteps))   cout << "50 % of the integration is performed" << endl;
        if(counter == floor(0.75*nSteps))   cout << "75 % of the integration is performed" << endl;
    }

    // Exploring the system using the Verlet integrator
    // Qualities of the system we will be exploring
    System solarSystemVerlet;
    solarSystemVerlet.addBody(0, 0, 0, 0, -M_PI*6e-6, 0, 1);
    solarSystemVerlet.addBody(1, 0, 0, 0, 2*M_PI, 0, 3e-6);

    Integrator verletsolver;
    //CelestialBody &twobodies = solarSystemVerlet.bodies[0];

    Printing printerv("Verlet");
    printerv.printingPosition(solarSystemVerlet);

    // Performing Verlet on the system
    for(int i = 0; i < nSteps; i++){
        verletsolver.Verlet(solarSystemVerlet, dt);
        printerv.printingPosition(solarSystemVerlet);
    }

    return 0;
}
