#include "KeplerClasses.h"


// ============================================== CLASS: CELESTIALBODY ============================================== //

// Initialising functions
CelestialBody::CelestialBody(){
}

//CelestialBody::CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass){
CelestialBody::CelestialBody(vec position, vec velocity, double mass){
    // Function defines a celestial body by position, velocity and mass from the input parameters

    this->position = position;
    this->velocity = velocity;
    this->mass     = mass;
}


// ============================================== CLASS: SYSTEM ===================================================== //

System::System(){
    nBodies = 0;
    dt = 1/12.0;
}


void System::addBody(vec position, vec velocity, double mass){
    CelestialBody body(position, velocity, mass);
    bodies[nBodies] = body;
    this->nBodies ++;
} // End addBody-function


void System::definingBodies(){
    // Creating system described by a vector defineSystem
    // with celestial bodies defined according to the CelestialBody class

    // The vector defineSystem is created to contain all position, velocity and mass information about the system
    defineSystem = zeros(7*nBodies);
    energyAngMom = zeros(4*nBodies);

    for(int i = 0; i < nBodies; i++){             // Defining the system by position, velocity and mass
        defineSystem[7*i]   = bodies[i].position[0];
        defineSystem[7*i+1] = bodies[i].position[1];
        defineSystem[7*i+2] = bodies[i].position[2];
        defineSystem[7*i+3] = bodies[i].velocity[0];
        defineSystem[7*i+4] = bodies[i].velocity[1];
        defineSystem[7*i+5] = bodies[i].velocity[2];
        defineSystem[7*i+6] = bodies[i].mass;
    }   // End for-loop defining the system
} // End definingBodies-function


void System::solver(double nYears){     // what is the function supposed to receive?? //

    double t = 0;

    ofstream posVelOut;
    ofstream enMom;
    posVelOut.open("posvel.dat");
    enMom.open("energyangmom.dat");

    while(t < nYears){
        RK4();

        // Writing position and velocity to file
        posVelOut << setiosflags(ios::showpoint | ios:: uppercase);
        posVelOut << setw(20) << setprecision(15) << defineSystem.t();
        // Writing energy and angular momentum to file
        enMom << setiosflags(ios::showpoint | ios:: uppercase);
        enMom << setw(20) << setprecision(15) << energyAngMom.t();

        t += dt;
    }
    posVelOut.close();
} // End solver-function


void System::RK4(){
    // Function performing the Runge-Kutta 4 method

    // Declaring slopes to perform RK4
    vec k0, k1, k2, k3, k4;

    k0 = zeros(7*nBodies);
    k1 = dt*diffEq(k0);
    k2 = dt*diffEq(k1/2.);
    k3 = dt*diffEq(k2/2.);
    k4 = dt*diffEq(k3);

    // Updating position and velocity after RK4
    defineSystem = defineSystem + (k1 + 2.*k2 + 2.*k3 + k4)/6.;

}   // Ending RK4-function

vec System::diffEq(vec k){
    // Function setting the differential equations

    // Creating temporary position and velocity vector, given by the "original" position and velocity,
    // possibly changed by a slope k found through the RK4-method
    vec posVel(7*nBodies);
    posVel = defineSystem + k;

    vec dotPosVel(7*nBodies);               // Creating derivative of position and velocity vector to return
    vec Fx = zeros(nBodies);                // Gravitational force in x-direction to compute vxdot
    vec Fy = zeros(nBodies);                // Gravitational force in y-direction to compute vydot
    vec Fz = zeros(nBodies);                // Gravitational force in z-direction to compute vzdot

    calculateForcesAndEnergy(posVel, Fx, Fy, Fz);

    // Finding the derivative of each velocity from the gravitational force
    for(int i = 0; i < nBodies; i++){
        dotPosVel[7*i+0] = posVel[7*i+3];               // setting xdot  = vx
        dotPosVel[7*i+1] = posVel[7*i+4];               // setting ydot  = vy
        dotPosVel[7*i+2] = posVel[7*i+5];               // setting zdot  = vz
        dotPosVel[7*i+3] = Fx[i]/defineSystem[7*i+6];   // setting vxdot = Fx/m
        dotPosVel[7*i+4] = Fy[i]/defineSystem[7*i+6];   // setting vydot = Fy/m
        dotPosVel[7*i+5] = Fz[i]/defineSystem[7*i+6];   // setting vzdot = Fz/m
    }   // Ending for-loop computing derivatives

    return dotPosVel;
}   // End of diffEq-function



/*
void System::calculateForcesAndEnergy(vec posVel, vec &Fx, vec &Fy, vec &Fz){
    // Function calculating forces and energy (and angular momentum!) for the system

    double G = 4*M_PI*M_PI;                 // Defining the gravitational constant in appropriate units
    Fx.zeros();
    Fy.zeros();
    Fz.zeros();

    for(int i = 0; i < nBodies; i++){
        for(int j=i+1; j < nBodies; j++){
            // Variables simplifying the calculations
            double dx     = posVel[7*i+0]-posVel[7*j+0];    // x-separation between celestial bodies
            double dy     = posVel[7*i+1]-posVel[7*j+1];    // y-separation between celestial bodies
            double dz     = posVel[7*i+2]-posVel[7*j+2];    // z-separation between celestial bodies
            double dr     = sqrt(dx*dx + dy*dy + dz*dz);    // separation between celestial bodies
            double factor = G*masses(j)*masses(i) / pow(dr,3);

            // Updating gravitational force and potential energy experienced by celestial object
            energyAngMom(1) += factor*dr;                   // Potential energy
            Fx[i] -= factor*dx;                             // Newton's law of gravity, x-direction
            Fy[i] -= factor*dy;                             // Newton's law of gravity, y-direction
            Fz[i] -= factor*dz;                             // Newton's law of gravity, z-direction
            Fx[j] += factor*dx;                             // Using N3 to find the force on the opposite object
            Fy[j] += factor*dy;                             // Using N3 to find the force on the opposite object
            Fz[j] += factor*dz;                             // Using N3 to find the force on the opposite object
        }   // Ending for-loop computing force and potential energy

        // Variables simplifying the calculations
        double pos = sqrt(pow(posVel[7*i+0],2) + pow(posVel[7*i+1],2) + pow(posVel[7*i+2],2)); // Position in CoM frame
        double v2  = pow(posVel[7*i+3],2) + pow(posVel[7*i+4],2) + pow(posVel[7*i+5],2);       // Velocity squared

        // Updating kinetic energy and angular momentum
        energyAngMom(3) += masses(i)*pos*sqrt(v2);                      // Angular momentum
        energyAngMom(0) += 0.5*masses(i)*v2;                            // Kinetic energy

    }   // Ending for-loop going over all celestial bodies

    energyAngMom(2) = sqrt(pow(energyAngMom(0),2) + pow(energyAngMom(1),2));   // Total energy

}   // Ending calculateForcesAndEnergy-function
*/
