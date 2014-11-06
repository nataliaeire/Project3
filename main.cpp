#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <system.h>
#include <integrator.h>
#include <printing.h>
#include <gaussiandeviate.h>
#include "time.h"

using namespace std;

// Functions possible to include in main function
void regularSystemRK4(double dt, double nSteps);
void regularSystemV(double dt, double nSteps);
void regularSystemVV(double dt, double nSteps);
void regularSystemVVadaptive(double nSteps);
void randomSystem(int numberOfObjects, double sphereRadius, double timeStep, double runningTime, bool smoothing);
void Mercury(double dt, double nSteps);

int main()
{
    // Integration specifications
    double dt       = 0.05;
    double T        = 500;
    double nSteps   = T/dt;
    int    numberOfObjects = 100;      // Number of celestial bodies for a random generation of a system
    double sphereRadius    = 20;       // ly

    // Running the code for special cases
    // Note that these functions were moved outside the main function simply to obtain a pretty, short main function,
    // as well as easily being able to run only some special cases, without having to comment away huge chunks
    // of the main function at a time
    //regularSystemRK4(dt, nSteps);       // Running the code using RK4
    //regularSystemV(dt, nSteps);         // Running the code using Verlet
    //regularSystemVV(dt, nSteps); // Running the code using Velocity Verlet
    //regularSystemVVadaptive(nSteps); // Running the code using Velocity Verlet
    // Mercury(dt, nSteps);             // Running the code for the GR case for Mercury
    randomSystem(numberOfObjects, sphereRadius, 0.005, 2, true);

    return 0;
}


void regularSystemRK4(double dt, double nSteps)
{ // Function analysing the given system using RK4
    // Qualities of the system we will be exploring are read from file
    fstream file("/home/ubu/FYS3150/projects/Project3/SunEarthNASA.txt",ios_base::in);
    if(!file.is_open()) {       // Printing error message about not being able to find file (at said location)
        cout << "Could not find file to open." << endl;
        exit(1);                                    // Cancelling the rest of the program
    }
    // Initialisation
    System      SolarSystem;                        // Preparing system
    Integrator  solving;                            // Preparing for allowing the system to develop
    Printing    printer("SunEarthNASARK4");         // Preparing for printing details about system to file

    SolarSystem.addSystem(file);                    // Creating system
    SolarSystem.conserveMomentum();                 // Ensuring momentum is conserved for the system
    printer.printingAll(SolarSystem);               // Printing intitial values to file

    //SolarSystem.sortBodiesIntoGroups();
    //SolarSystem.calculateForcesAdaptively(8);
    // NNNNNNNNNNNNNNNNNNNNNBBBBBBBBBBBBBBBBBBBBBBB!

    int printNFrames = 1;                           // Counter for printing only each n'th frame

    double start = clock();
    // Performing RK4 on the system
    for(int i = 0; i < nSteps; i++){
        solving.RK4(SolarSystem, dt);               // Solving the problem using the RK4-method
        printer.printingAll(SolarSystem, i, printNFrames);    // Printing everything to file

        // Printing a message to screen to let the user know how far the program has come
        if(i % int(1e3) == 0)     cout << 100*((double)i) / nSteps << " % of the RK4 integration is performed" << endl;
    }
    double finish = clock();
    double operationTime = (finish - start)/(double) CLOCKS_PER_SEC/nSteps; // Calculating time in seconds
    cout << "Operation time pr. RK4 time step: " << operationTime << " s" << endl << endl;

} // End of regularSystemRK4-function


void regularSystemV(double dt, double nSteps)
{ // Function analysing the given system using Verlet
    // Qualities of the system we will be exploring are read from file
    fstream file("/home/ubu/FYS3150/projects/Project3/SunEarthNASA.txt",ios_base::in);
    if(!file.is_open()) {       // Printing error message about not being able to find file (at said location)
        cout << "Could not find file to open." << endl;
        exit(1);                                    // Cancelling the rest of the program
    }

    // Initialisation
    System      solarSystemVerlet;
    Integrator  verletsolver;
    Printing    printerv("SunEarthNASAV");

    solarSystemVerlet.addSystem(file);
    solarSystemVerlet.conserveMomentum();
    printerv.printingAll(solarSystemVerlet);

    int printNFrames = 1;                           // Counter for printing only each n'th frame

    double start = clock();
    // Performing Verlet on the system
    for(int i = 0; i < nSteps; i++){
        verletsolver.Verlet(solarSystemVerlet, dt);
        printerv.printingAll(solarSystemVerlet, i, printNFrames);

        // Printing a message to screen to let the user know how far the program has come
        if(i % int(1e3) == 0)     cout << 100*((double)i) / nSteps << " % of the Verlet integration is performed" << endl;
    }
    double finish = clock();
    double operationTime = (finish - start)/(double) CLOCKS_PER_SEC/nSteps; // Calculating time in seconds
    cout << "Operation time pr. Verlet time step: " << operationTime << " s" << endl << endl;

} // End of regularSystemV-function

void regularSystemVV(double dt, double nSteps)
{ // Function analysing the given system using Velocity Verlet
    // Qualities of the system we will be exploring are read from file
    fstream file("/home/ubu/FYS3150/projects/Project3/SunEarthNASA.txt",ios_base::in);
    if(!file.is_open()) {       // Printing error message about not being able to find file (at said location)
        cout << "CCould not find file to open." << endl;
        exit(1);                                    // Cancelling the rest of the program
    }

    // Initialisation
    System      solarSystemVelocityVerlet;
    Integrator  velocityverletsolver;
    Printing    printervv("solarsystemNASAVV");

    solarSystemVelocityVerlet.addSystem(file);
    solarSystemVelocityVerlet.conserveMomentum();
    printervv.printingAll(solarSystemVelocityVerlet);

    int printNFrames = 1;                           // Counter for printing only each n'th frame

    double start = clock();
    // Performing Verlet on the system
    for(int i = 0; i < nSteps; i++){
        velocityverletsolver.VelocityVerlet(solarSystemVelocityVerlet, dt);
        printervv.printingAll(solarSystemVelocityVerlet, i, printNFrames);

        // Printing a message to screen to let the user know how far the program has come
        if(i % int(1e3) == 0)     cout << 100*((double)i) / nSteps << " % of the Velocity Verlet integration is performed" << endl;
    }
    double finish = clock();
    double operationTime = (finish - start)/(double) CLOCKS_PER_SEC/nSteps; // Calculating time in seconds
    cout << "Operation time pr. Velocity Verlet time step: " << operationTime << " s" << endl << endl;

} // End of regularSystemVV-function


void regularSystemVVadaptive( double nSteps)
{ // Function analysing the given system using Velocity Verlet
    // Qualities of the system we will be exploring are read from file
    fstream file("/home/ubu/FYS3150/projects/Project3/solarsystemNASA.txt",ios_base::in);
    if(!file.is_open()) {       // Printing error message about not being able to find file (at said location)
        cout << "CCould not find file to open." << endl;
        exit(1);                                    // Cancelling the rest of the program
    }

    // Initialisation
    System      solarSystemVelocityVerlet;
    Integrator  velocityverletsolver;
    Printing    printervv("solarsystemNASAVV");

    solarSystemVelocityVerlet.addSystem(file);
    solarSystemVelocityVerlet.conserveMomentum();
    printervv.printingAll(solarSystemVelocityVerlet);

    int printNFrames = 1;                           // Counter for printing only each n'th frame

    double start = clock();
    // Performing Verlet on the system
    for(int i = 0; i < nSteps; i++){
        velocityverletsolver.adaptiveVelocityVerlet(solarSystemVelocityVerlet);
        printervv.printingAll(solarSystemVelocityVerlet, i, printNFrames);

        // Printing a message to screen to let the user know how far the program has come
        if(i % int(1e3) == 0)     cout << 100*((double)i) / nSteps << " % of the Velocity Verlet integration is performed" << endl;
    }
    double finish = clock();
    double operationTime = (finish - start)/(double) CLOCKS_PER_SEC/nSteps; // Calculating time in seconds
    cout << "Operation time pr. Velocity Verlet time step: " << operationTime << " s" << endl << endl;

} // End of regularSystemVV-function


void Mercury(double dt, double nSteps)
{ // Function performing the evolution of Mercury using a GR contribution
    // Qualities of the system we will be exploring are read from file
    fstream MercuryFile("/home/ubu/FYS3150/projects/Project3/Mercury.txt",ios_base::in);
    if(!MercuryFile.is_open()) {    // Printing error message about not being able to find file (at said location)
        cout << "Could not find file to open." << endl;
        exit(1);                                    // Cancelling the rest of the program
    } // End if-statement

    // Intialisation
    System      MercurySystem;                      // Preparing system
    Integrator  solvingMercury;                     // Preparing for allowing the system to develop
    Printing    printerMercury("Mercury");          // Preparing for printing details about system to file

    MercurySystem.addSystem(MercuryFile);           // Creating system
    MercurySystem.conserveMomentum();               // Ensuring momentum is conserved for the system

    vec3 old(1.,1.,1.);
    vec3 superold(0.5,0.5,0.5);

    // Performing RK4 on the system
    for(int i = 0; i < nSteps; i++){
        solvingMercury.RK4GR(MercurySystem, dt);    // Solving the problem using the RK4-method
        CelestialBody Mercury = MercurySystem.bodies[1];

        // Print previous position to file, if the system was at perihelion at the previous step
        if( old.length() < Mercury.position.length() && old.length() < superold.length() ){
            printerMercury.printing3Vector(old,"position");
        } // End if-statement

        superold = old;
        old = Mercury.position;

        // Printing a message to screen to let the user know how far the program has come
        if(i % int(1e5) == 0)     cout << 100*((double)i) / nSteps << " % of the integration is performed" << endl;

    } // Ending for-loop
} // End of Mercury-function


void randomSystem(int numberOfObjects, double sphereRadius, double timeStep, double runningTime, bool smoothing)
{ // Doing calculations for a randomly generated system
    // Initialisation
    System      system;
    Integrator  solvingSystem;
    Printing    printingSystem("RandomSystem100bodies");

    system.smoothing = smoothing;
    system.addRandomSystem(numberOfObjects,sphereRadius);
    printingSystem.printingPositionXYZ(system);

    double  time = 0;
    double  nextPrintTime = 0;
    int     counter = 1;    // Counter to keep track of time steps in output file
    int     numberOfTimestepsComputed = 0;

    while(time < runningTime){
        numberOfTimestepsComputed++;
        solvingSystem.adaptiveVelocityVerlet(system);
        printingSystem.printingEnergyAngMom(system,true);

        time += 8.*solvingSystem.adaptiveDt();
        //time += timeStep;

        if(time > nextPrintTime){
            printingSystem.printingPositionXYZ(system, counter);
            nextPrintTime += 0.02*runningTime;
            cout << 100*(time/runningTime) << " % of the Velocity Verlet integration is performed, currently at " << numberOfTimestepsComputed << " timesteps, with " << solvingSystem.adaptiveDt() <<endl;
            counter++;
        } // Ending if-statement
    } // End while-loop
} // End of randomSystem-function
