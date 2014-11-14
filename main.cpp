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
void randomSystemNonAdaptive(int numberOfObjects, double sphereRadius, double timeStep, double runningTime, bool smoothing);
void randomSystemAdaptive(int numberOfObjects, double sphereRadius, double runningTime, bool smoothing);
void Mercury(double dt, double nSteps);

int main()
{
    // Integration specifications for Solar System

    double dt       = 1e-5;
    double T        = 1;
    double nSteps   = T/dt;
    //double nSteps   = (T+dt)/dt; // Something like this fixes the problem with the Earth not completing its orbit in one year

    // Specifics for cold collapse
    double timeStep         = 0.005;    // tcrunch
    double runningTime      = 5;        // tcrunch
    double sphereRadius     = 20;       // ly
    int    numberOfObjects  = 750;      // Number of celestial bodies for a random generation of a system
    bool   smoothing        = true;

    // Running the code for special cases
    // Note that these functions were moved outside the main function simply to obtain a pretty, short main function,
    // as well as easily being able to run only some special cases, without having to comment away huge chunks
    // of the main function at a time

    regularSystemRK4(dt, nSteps);       // Running the code using RK4
    regularSystemV(dt, nSteps);         // Running the code using Verlet
    //regularSystemVV(dt, nSteps);        // Running the code using Velocity Verlet
    //regularSystemVVadaptive(nSteps);    // Running the code using Velocity Verlet
    //Mercury(dt, nSteps);             // Running the code for the GR case for Mercury
    //randomSystemNonAdaptive(numberOfObjects, sphereRadius, timeStep, runningTime, smoothing);
    //randomSystemAdaptive(numberOfObjects, sphereRadius, runningTime, smoothing);

    return 0;
}


void regularSystemRK4(double dt, double nSteps)
{ // Function analysing the given system using RK4
    // Qualities of the system we will be exploring are read from file
    fstream file("/uio/hume/student-u67/kineoha/FYS3150/Project3/SunEarthNASA.txt",ios_base::in);
    if(!file.is_open()) {       // Printing error message about not being able to find file (at said location)
        cout << "Could not find file to open." << endl;
        exit(1);                                    // Cancelling the rest of the program
    }
    // Initialisation
    System      SolarSystem;                        // Preparing system
    Integrator  solving;                            // Preparing for allowing the system to develop
<<<<<<< HEAD

    Printing    printer("SolarSystemNASARK4dtexp_minus6pnf1e4");         // Preparing for printing details about system to file
=======
    Printing    printer("SolarSystemNASAbubblebath");         // Preparing for printing details about system to file
>>>>>>> 3444baf9ab1d8f04243cd872bf13b10afe935499

    //SolarSystem.addSystem(file);                    // Creating system
    SolarSystem.addBody(0,0,0,0,0,0,1);
    SolarSystem.addBody(1,0,0, 0, 2*M_PI, 0, 3e-6);
    SolarSystem.conserveMomentum();                 // Ensuring momentum is conserved for the system
    printer.printingAll(SolarSystem);               // Printing intitial values to file

    int printNFrames = 1e3;                           // Counter for printing only each n'th frame


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
    fstream file("/uio/hume/student-u67/kineoha/FYS3150/Project3/SunEarthNASA.txt",ios_base::in);
    if(!file.is_open()) {       // Printing error message about not being able to find file (at said location)
        cout << "Could not find file to open." << endl;
        exit(1);                                    // Cancelling the rest of the program
    }

    // Initialisation
    System      solarSystemVerlet;
    Integrator  verletsolver;
<<<<<<< HEAD
    Printing    printerv("SunEarthNASAVerletoneyeardtexp_minus1pnfexp0");
=======
    Printing    printerv("SunEarthNASAVerletbubblebath");
>>>>>>> abec4eda4e6eabc9c5adbf238a2f390a111aeae7


    //solarSystemVerlet.addSystem(file);
    solarSystemVerlet.addBody(0,0,0,0,0,0,1);
    solarSystemVerlet.addBody(1,0,0, 0, 2*M_PI, 0, 3e-6);
    solarSystemVerlet.conserveMomentum();
    printerv.printingAll(solarSystemVerlet);

    int printNFrames = 1e3;                          // Counter for printing only each n'th frame

    double start = clock();
    // Performing Verlet on the system
    for(int i = 0; i < nSteps; i++){
        verletsolver.Verlet(solarSystemVerlet, dt);
        printerv.printingAll(solarSystemVerlet, i, printNFrames);

        // Printing a message to screen to let the user know how far the program has come
        if(i % int(1e3) == 0)     cout << 100*((double)i) / nSteps << " % of the Verlet integration is performed\r";
    }
    double finish = clock();
    double operationTime = (finish - start)/(double) CLOCKS_PER_SEC/nSteps; // Calculating time in seconds
    cout << "Operation time pr. Verlet time step: " << operationTime << " s" << endl << endl;

} // End of regularSystemV-function

void regularSystemVV(double dt, double nSteps)
{ // Function analysing the given system using Velocity Verlet
    // Qualities of the system we will be exploring are read from file
    fstream file("/uio/hume/student-u67/kineoha/FYS3150/Project3/SunEarthNASA.txt",ios_base::in);
    if(!file.is_open()) {       // Printing error message about not being able to find file (at said location)
        cout << "Could not find file to open." << endl;
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
    fstream file("/uio/hume/student-u67/kineoha/FYS3150/Project3/solarsystemNASA.txt",ios_base::in);
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
    fstream MercuryFile("/uio/hume/student-u67/kineoha/FYS3150/Project3/Mercury.txt",ios_base::in);
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


void randomSystemNonAdaptive(int numberOfObjects, double sphereRadius, double timeStep, double runningTime, bool smoothing)
{ // Doing calculations for a randomly generated system
    // Initialisation
    System      system;
    Integrator  solvingSystem(1);
    Printing    printingSystem("test");

    system.smoothing = smoothing;
    system.addRandomSystem(numberOfObjects,sphereRadius);
    printingSystem.printingPositionXYZ(system);

    double  time = 0;
    double  nextPrintTime = 0;
    int     counter = 1;                                // Counter to keep track of time steps in output file
    int     numberOfTimestepsComputed = 0;

    double start = clock();
    // Evolving system
    while(time < runningTime){
        numberOfTimestepsComputed++;
        solvingSystem.VelocityVerlet(system, timeStep);

        time += timeStep;

        if(time > nextPrintTime){
            printingSystem.printingPositionXYZ(system, counter);
            printingSystem.printingEnergyAngMom(system,true);
            nextPrintTime += 0.002*runningTime;
            cout << 100*(time/runningTime) << " % of the Velocity Verlet integration is performed, currently at "
                 << numberOfTimestepsComputed << " timesteps, with time step " << solvingSystem.adaptiveDt()
                 << "." << endl;
            counter++;
        } // Ending if-statement
    } // End while-loop
    double finish = clock();
    double operationTime = (finish - start)/(double) CLOCKS_PER_SEC; // Calculating time in seconds
    cout << "Operation time: " << operationTime << " s" << endl << endl;
    cout << "Parallelisation:" << solvingSystem.numThreads << endl;

} // End of randomSystem-function


void randomSystemAdaptive(int numberOfObjects, double sphereRadius, double runningTime, bool smoothing)
{ // Doing calculations for a randomly generated system
    // Initialisation
    System      system;
    Integrator  solvingSystem(4);
    Printing    printingSystem("VirialTest");

    // Variables to be used in the integration
    double  time = 0;
    double  nextPrintTime = 0;
    int     counter = 1;                                // Counter to keep track of time steps in output file
    int     numberOfTimestepsComputed = 0;

    // Creating a specific system and saving intials to file
    system.smoothing = smoothing;
    system.addRandomSystem(numberOfObjects,sphereRadius);
    printingSystem.printingPosition(system, true);
    printingSystem.printingPositionXYZ(system, counter);

    double start = clock();
    // Evolving system
    while(time < runningTime){
        numberOfTimestepsComputed++;
        solvingSystem.adaptiveVelocityVerlet(system);

        time += 8.*solvingSystem.adaptiveDt();

        if(time > nextPrintTime){
            printingSystem.printingAll(system, counter, true);
            nextPrintTime += 0.002*runningTime;
            cout << 100*(time/runningTime) << " % of the Velocity Verlet integration is performed, currently at "
                 << numberOfTimestepsComputed << " timesteps, with time step " << solvingSystem.adaptiveDt()
                 << "." << endl;
            counter++;
        } // Ending if-statement
    } // End while-loop
    double finish = clock();
    double operationTime = (finish - start)/(double) CLOCKS_PER_SEC; // Calculating time in seconds
    cout << endl << "Operation time: " << operationTime << " s" << endl;
    cout << "Parallelisation/Number of threads: " << solvingSystem.numThreads << endl;

} // End of randomSystem-function
