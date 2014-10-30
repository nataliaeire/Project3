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
void findError();
void regularSystemRK4(double dt, double nSteps);
void regularSystemV(double dt, double nSteps);
void regularSystemVV(double dt, double nSteps);
void regularSystemVVadaptive(double nSteps);
void randomSystem(int numberOfObjects, double sphereRadius, double timeStep, double runningTime, bool smoothing);
void Mercury(double dt, double nSteps);

int main()
{
    // findError() finds the error introduced by increasing the time step.
    // Should be commented away if not being used, as the rest of the main function is ignored if it runs succesfully.
    // findError();

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
    randomSystem(250, sphereRadius, 0.01, 5, true);

    /*
    System system;
    Printing printing("initialpositions");
    system.addRandomSystem(1000,20);
    printing.printingPositionXYZ(system);
    */

    return 0;
}


void findError()
{ // Function aiming to find a quantitative error for different time steps
    // Creating file to write the errors to
    ofstream errorFile;
    errorFile.open("Error.txt");

    // Preparing system
    fstream     changingFile("/home/ubu/FYS3150/projects/Project3/SunEarthNASA.txt",ios_base::in);

    if(!changingFile.is_open()) {
        cout << "Could not find file: " << endl;
        exit(1);
    }

    System      initialSystem;
    initialSystem.addSystem(changingFile);
    initialSystem.conserveMomentum();

    // Quantities needed in the calculation
    arma::vec   errorRK4;
    arma::vec   errorVV;
    errorRK4.zeros(7);
    errorVV.zeros(7);
    vec3        errorVecRK4;
    vec3        errorVecVV;
    vec3        rTrueRK4;
    vec3        rTrueVV;
    vec3        rDtRK4;
    vec3        rDtVV;
    double      dtChanging = 1e-7;
    double      changingNSteps;

    initialSystem.calculateForcesAndEnergy();
    double      initialEnergy = initialSystem.totalEnergy();

    // Computing error for different time step values
    for(int i = 0; i<7; i++){
        dtChanging     *= 10;
        changingNSteps  = 1./dtChanging;

        System changingSystemRK4 = initialSystem;
        System changingSystemVV  = initialSystem;
        Integrator  changingSolverRK4;
        Integrator  changingSolverVV;

        // Evolving the system for a year
        for(int j = 0; j < changingNSteps; j++){
            changingSolverRK4.RK4(changingSystemRK4, dtChanging);
            changingSolverVV.Verlet(changingSystemVV,dtChanging);
        } // Ending for-loop

        CelestialBody earthRK4 = changingSystemRK4.bodies[1];
        CelestialBody earthVV  = changingSystemVV.bodies[1];
        rDtRK4 = earthRK4.position;
        rDtVV  = earthVV.position;

        // Defining the true value of the position of the Earth to be one where the time step was small
        if(i==0){
            rTrueRK4 = rDtRK4;
            rTrueVV  = rDtVV;
        } // Ending if-statement

        // Computing the error
        errorVecRK4 = rTrueRK4 - rDtRK4;
        errorVecVV  = rTrueVV  - rDtVV;
        //errorRK4(i) = errorVecRK4.length();
        //errorVV(i)  = errorVecVV.length();
        errorRK4(i) = changingSystemRK4.totalEnergy() - initialEnergy;
        errorVV(i)  = changingSystemVV.totalEnergy()  - initialEnergy;

    } // Ending for-loop

    // Printing error tofile
    errorFile << errorRK4.t();
    errorFile << errorVV.t() << endl;

    errorFile.close();
    exit(0);            // Shut down the rest of the main function
} // End of findError-function


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
    SolarSystem.calculateForcesAdaptively(8);

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
        velocityverletsolver.adaptiveVelocityVerlet(solarSystemVelocityVerlet, i);
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
{
    // Initialisation
    System      system;
    Integrator  solvingSystem;
    Printing    printingSystem("RandomSystem");

    system.addRandomSystem(numberOfObjects,sphereRadius);
    system.smoothing = smoothing;
    printingSystem.printingPositionXYZ(system);

    double time = 0;
    double nextPrintTime = 0;

    while(time < runningTime){
        solvingSystem.adaptiveVelocityVerlet(system,timeStep,true);


        time += solvingSystem.adaptiveDt();

        if(time > nextPrintTime){
            printingSystem.printingPositionXYZ(system);
            nextPrintTime += runningTime/500;
        } // Ending if-statement

        /*
        CelestialBody body1 = system.bodies[0];
        CelestialBody body2 = system.bodies[1];
        cout << body1.mass << "     " << body2.mass << endl;
        cout << body1.force << endl;
        cout << body2.force << endl << endl;
        */
        if((i+1) % 8 == 0)      printingSystem.printingEnergyAngMom(system);

        // Printing a message to screen to let the user know how far the program has come
        if(i % int(10) == 0)     cout << 100*((double)i) / n << " % of the Velocity Verlet integration is performed" << endl;
    }
} // End of randomSystem-function
