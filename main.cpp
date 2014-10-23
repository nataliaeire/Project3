#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <system.h>
#include <integrator.h>
#include <printing.h>

using namespace std;

void findError();
void regularSystemRK4(double dt, double nSteps);
void regularSystemVV(double dt, double nSteps);
void Mercury(double dt, double nSteps);

int main()
{
    // findError() finds the error introduced by increasing the time step.
    // Should be commented away if not being used, as the rest of the main function is ignored if it runs succesfully.
    // findError();

    // Integration specifications
    double dt       = 1e-6;
    double T        = 5;
    double nSteps   = T/dt;

    // Running the code for special cases
    // Note that these functions were moved outside the main function simply to obtain a pretty, short main function,
    // as well as easily being able to run only some special cases, without having to comment away huge chunks
    // of the main function at a time
    regularSystemRK4(dt, nSteps);       // Running the code using RK4
    regularSystemVV(dt, nSteps);        // Running the code using VV
    // Mercury(dt, nSteps);             // Running the code for the GR case for Mercury

    return 0;
}


void findError()
{
    // Creating file to write the errors to
    ofstream errorFile;
    errorFile.open("Error.txt");

    // Preparing system
    fstream     changingFile("/uio/hume/student-u81/natales/Project3/Project3/SunEarthNASA.txt",ios_base::in);
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

    // Computing error for different time step values
    for(int i = 0; i<7; i++){
        dtChanging     *= 10;
        changingNSteps  = 1/dtChanging;

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
        errorRK4(i) = errorVecRK4.length();
        errorVV(i)  = errorVecVV.length();
    } // Ending for-loop

    // Printing error tofile
    errorFile << errorRK4.t();
    errorFile << errorVV.t() << endl;

    errorFile.close();

    exit(0);        // Shut down the rest of the main function
} // End of findError-function


void regularSystemRK4(double dt, double nSteps)
{
    // Qualities of the system we will be exploring are read from file
    fstream file("/uio/hume/student-u81/natales/Project3/Project3/SunEarthNASA.txt",ios_base::in);

    // Initialisation
    System      SolarSystem;                 // Preparing system
    Integrator  solving;                     // Preparing for allowing the system to develop
    Printing    printer("SolarSystemRK4");   // Preparing for printing details about system to file

    SolarSystem.addSystem(file);             // Creating system
    SolarSystem.conserveMomentum();          // Ensuring momentum is conserved for the system
    printer.printingAll(SolarSystem);        // Printing intitial values to file

    int printNFrames = 1e3;                  // Counter for printing only each n'th frame

    // Performing RK4 on the system
    for(int i = 0; i < nSteps; i++){
        solving.RK4(SolarSystem, dt);        // Solving the problem using the RK4-method
        printer.printingAll(SolarSystem, i, printNFrames);    // Printing everything to file

        // Printing a message to screen to let the user know how far the program has come
        if(i % 100000 == 0)     cout << 100*((double)i) / nSteps << " % of the RK4 integration is performed" << endl;
    }
} // End of regularSystemRK4-function


void regularSystemVV(double dt, double nSteps)
{
    // Qualities of the system we will be exploring are read from file
    fstream file("/uio/hume/student-u81/natales/Project3/Project3/SunEarthNASA.txt",ios_base::in);

    // Initialisation
    System      solarSystemVerlet;
    Integrator  verletsolver;
    Printing    printerv("SolarSystemVV");

    solarSystemVerlet.addSystem(file);
    solarSystemVerlet.conserveMomentum();
    printerv.printingAll(solarSystemVerlet);

    int printNFrames = 1e3;                 // Counter for printing only each n'th frame

    // Performing Verlet on the system
    for(int i = 0; i < nSteps; i++){
        verletsolver.Verlet(solarSystemVerlet, dt);
        printerv.printingAll(solarSystemVerlet, i, printNFrames);

        // Printing a message to screen to let the user know how far the program has come
        if(i % 100000 == 0)     cout << 100*((double)i) / nSteps << " % of the Verlet integration is performed" << endl;
    }
} // End of regularSystemVV-function


void Mercury(double dt, double nSteps)
{
    // Qualities of the system we will be exploring are read from file
    fstream MercuryFile("/home/ubu/FYS3150/projects/Project3/MercuryInitials.txt",ios_base::in);

    // Intialisation
    System      MercurySystem;                      // Preparing system
    Integrator  solvingMercury;                     // Preparing for allowing the system to develop
    Printing    printerMercury("Mercury");          // Preparing for printing details about system to file

    MercurySystem.addSystem(MercuryFile);           // Creating system
    MercurySystem.conserveMomentum();               // Ensuring momentum is conserved for the system

    vec3 old(1.,1.,1.);
    vec3 superold(0.5,0.5,0.5);
    //vec3 oldvelocity(0.,0.,0.);

    // Performing RK4 on the system
    for(int i = 0; i < nSteps; i++){
        solvingMercury.RK4GR(MercurySystem, dt);    // Solving the problem using the RK4-method
        CelestialBody Mercury = MercurySystem.bodies[1];

        // Print previous position to file, if the system was at perihelion at the previous step
        if( old.length() < Mercury.position.length() && old.length() < superold.length() ){
            printerMercury.printing3Vector(old,"position");
            //cout << oldvelocity << endl;
        } // End if-statement

        superold = old;
        old = Mercury.position;
        //oldvelocity = Mercury.velocity;

        // Printing a message to screen to let the user know how far the program has come
        if(i % 10000 == 0)     cout << 100*((double)i) / nSteps << " % of the integration is performed" << endl;

    } // Ending for-loop
} // End of Mercury-function
