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
    double dt       = 0.000001;
    double T        = 10;
    double nSteps   = T/dt;


    /*
    // Qualities of the system we will be exploring are read from file
    // Note that file directory has to be changed accordingly for every computer
    fstream file("/uio/hume/student-u81/natales/Project3/Project3/hei.txt",ios_base::in);
    fstream file2("/uio/hume/student-u81/natales/Project3/Project3/hei.txt",ios_base::in);

    // Intialisation
    System      SunEarthSystem;                 // Preparing system
    Integrator  solving;                        // Preparing for allowing the system to develop
    Printing    printer("SunEarth");            // Preparing for printing details about system to file

    SunEarthSystem.addSystem(file);             // Creating system
    SunEarthSystem.conserveMomentum();          // Ensuring momentum is conserved for the system
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
    Integrator verletsolver;
    Printing printerv("Verlet");

    solarSystemVerlet.addSystem(file2);
    solarSystemVerlet.conserveMomentum();
    printerv.printingAll(solarSystemVerlet);

    // Performing Verlet on the system
    for(int i = 0; i < nSteps; i++){
        verletsolver.Verlet(solarSystemVerlet, dt);
        printerv.printingAll(solarSystemVerlet);
    }
    */


    // ============================== MERCURY ============================== //
    fstream MercuryFile("/uio/hume/student-u81/natales/Project3/Project3/MercuryInitials.txt",ios_base::in);

    // Intialisation
    System      MercurySystem;                      // Preparing system
    Integrator  solvingMercury;                     // Preparing for allowing the system to develop
    Printing    printerMercury("Mercury");          // Preparing for printing details about system to file

    MercurySystem.addSystem(MercuryFile);           // Creating system
    MercurySystem.conserveMomentum();               // Ensuring momentum is conserved for the system
    printerMercury.printingAll(MercurySystem);      // Printing intitial values to file

    int counterM = 1;                               // Counting parameter to print message to screen inside for-loop

    // Performing RK4 on the system
    for(int i = 0; i < nSteps; i++){
        solvingMercury.RK4GR(MercurySystem, dt);    // Solving the problem using the RK4-method
        printerMercury.printingAll(MercurySystem);  // Printing everything to file
        counterM ++;

        // Printing a message to screen to let the user know how far the program has come
        if(counterM == floor(0.25*nSteps))   cout << "25 % of the integration is performed" << endl;
        if(counterM == floor(0.50*nSteps))   cout << "50 % of the integration is performed" << endl;
        if(counterM == floor(0.75*nSteps))   cout << "75 % of the integration is performed" << endl;
    }

    return 0;
}
