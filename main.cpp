#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <system.h>
#include <integrator.h>
#include <printing.h>

using namespace std;

//TODO
int main()
{
    // Integration specifications
    double dt       = 1e-6;
    double T        = 100;
    double nSteps   = T/dt;



    // Qualities of the system we will be exploring are read from file
    // Note that file directory has to be changed accordingly for every computer
    fstream file("/home/ubu/FYS3150/projects/Project3/SunEarthNASA.txt",ios_base::in);
    fstream file2("/home/ubu/FYS3150/projects/Project3/SunEarthNASA.txt",ios_base::in);

    // Intialisation
    System      SolarSystem;                 // Preparing system
    Integrator  solving;                     // Preparing for allowing the system to develop
    Printing    printer("SunEarthNASA100years100pointsperyear");      // Preparing for printing details about system to file

    SolarSystem.addSystem(file);             // Creating system
    SolarSystem.conserveMomentum();          // Ensuring momentum is conserved for the system
    printer.printingAll(SolarSystem);        // Printing intitial values to file

    int print_n_frames = 1e4;                // Counter for printing only each n'th frame

    // Performing RK4 on the system
    for(int i = 0; i < nSteps; i++){
        solving.RK4(SolarSystem, dt);        // Solving the problem using the RK4-method
        printer.printingAll(SolarSystem, i, print_n_frames);    // Printing everything to file

        // Printing a message to screen to let the user know how far the program has come
        if(i % 1000 == 0)     cout << 100*((double)i) / nSteps << " % of the RK4 integration is performed" << endl;
    }

    // Exploring the system using the Verlet integrator
    // Qualities of the system we will be exploring

    System solarSystemVerlet;
    Integrator verletsolver;
    Printing printerv("VerletSunEarthNASA100years100pointsperyear");

    solarSystemVerlet.addSystem(file2);
    solarSystemVerlet.conserveMomentum();
    printerv.printingAll(solarSystemVerlet);

    // Performing Verlet on the system
    for(int i = 0; i < nSteps; i++){
        verletsolver.Verlet(solarSystemVerlet, dt);
        printerv.printingAll(solarSystemVerlet, i, print_n_frames);

        // Printing a message to screen to let the user know how far the program has come
        if(i % 1000 == 0)     cout << 100*((double)i) / nSteps << " % of the Verlet integration is performed" << endl;
    }
/*


    // ============================== MERCURY ============================== //
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

*/
    return 0;
}
