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

    // Qualities of the system we will be exploring
    System SunEarthSystem;
    SunEarthSystem.addBody(0, 0, 0, 0, -M_PI*6e-6, 0, 1);
    SunEarthSystem.addBody(1, 0, 0, 0, 2*M_PI, 0, 3e-6);

    Integrator solving;

    Printing printer("SunEarth");
    printer.printingAll(SunEarthSystem);        // Printing intitial values to file

    int counter = 1;                            // Counting parameter to print message to screen inside for-loop

    fstream file("filename.dat",ios_base::in);  // THIS FUNCTION NEEDS TO ACTUALLY BE WRITTEN //

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

    return 0;
}
