#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace arma;


// ============================================== CLASS: CELESTIALBODY ============================================== //
class CelestialBody{
public:
    char*  name;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double mass;

    CelestialBody();
    CelestialBody(double x, double y, double vx, double vy, double mass);
};  // End of CelestialBody class declaration


// Initialising functions
CelestialBody::CelestialBody(){
}

CelestialBody::CelestialBody(double x, double y, double vx, double vy, double mass){
    // Function defines a celestial body by position, velocity and mass from the input parameters

    //this->name = *name;
    this->x    = x;
    this->y    = y;
    //this->z    = z;
    this->vx   = vx;
    this->vy   = vy;
    //this->vz   = vz;
    this->mass = mass;
}


// ============================================== CLASS: SYSTEM ===================================================== //
class System {
public:
    int nBodies;
    CelestialBody* bodies;

    // Initialization
    System(int n);
    vec defineSystem;

    // Initialization functions
    void twoBodies();
    void threeBodies();
    void definingBodies();

    // Solving functions
    void RK4(double h, ofstream* posvel, ofstream* enmom);
    vec  diffEq(vec, ofstream*);

};  // End of System class declariation

System::System(int n){
    // Function declaring a system with a number (nBodies) of celestial bodies, each of which should be
    // defined by the way celestial bodies are declared in the CelestialBody class

    this->nBodies = n;
    this->bodies  = new CelestialBody[nBodies];

    for (int i = 0; i < nBodies; i ++){
        this->bodies[i] = CelestialBody();
    }
}

void System::twoBodies(){
    // Function defining a two-body system containing the Sun and the Earth

    // Sun
    //this->bodies[0].name = "Sun";
    this->bodies[0].x    = 0.0;
    this->bodies[0].y    = 0.0;
    //this->bodies[0].z    = 0.0;
    this->bodies[0].vx   = 0.0;
    this->bodies[0].vy   = -M_PI*6e-6;      // Initial velocity to allow the total momentum of the system to be zero
    //this->bodies[0].vz   = 0.0;
    this->bodies[0].mass = 1.0;

    // Earth
    //this->bodies[1].name = "Earth";
    this->bodies[1].x    = 1.0;
    this->bodies[1].y    = 0.0;
    this->bodies[1].vx   = 2*M_PI;
    this->bodies[1].vy   = 0.0;
    this->bodies[1].mass = 3e-6;
}

void System::threeBodies(){
    // Function defining a three-body system containing the Sun, the Earth and Jupiter

    // Sun
    this->bodies[0].mass = 1.0;
    this->bodies[0].x    = 0.0;
    this->bodies[0].y    = 0.0;
    this->bodies[0].vx   = 0.0;
    this->bodies[0].vy   = -M_PI*6e-6;      // Initial velocity to allow the total momentum of the system to be zero

    // Earth
    this->bodies[1].mass = 3e-6;
    this->bodies[1].x    = 1.0;
    this->bodies[1].y    = 0.0;
    this->bodies[1].vx   = 0.0;
    this->bodies[1].vy   = 2*M_PI;

    // Jupiter
    this->bodies[2].mass = 1e-3;
    this->bodies[2].x    = 0.0;
    this->bodies[2].y    = 5.20;
    this->bodies[2].vx   = 0.0;
    this->bodies[2].vy   = 0.0;
}

void System::definingBodies(){
    // Creating system described by a vector defineSystem
    // with celestial bodies defined according to the CelestialBody class

    // The vector defineSystem is created to contain all position, velocity and mass information about the system
    defineSystem(5*2); //this should be 5*nBodies

    for(int i = 0; i < 2; i++){             // Defining the system by position, velocity and mass
        defineSystem[5*i]   = bodies[i].x;
        defineSystem[5*i+1] = bodies[i].y;
        defineSystem[5*i+2] = bodies[i].vx;
        defineSystem[5*i+3] = bodies[i].vy;
        defineSystem[5*i+4] = bodies[i].mass;
    }   // End for-loop defining the system
}


/*
// This is what a function declaring a system with n bodies should look like.
void System::nBodies(int n){
    for (int i=1..n)
        this->bodies[i].x = rand();
}*/

/*
vec System::diffEq(vec k, ofstream* enmom){
    // Function setting the differential equations

    // Creating temporary position and velocity vector, given by the "original" position and velocity,
    // possibly changed by a slope k found through the RK4-method
    vec posVel(5*nBodies);
    posVel = defineSystem + k;

    vec dotPosVel(5*nBodies);               // Creating derivative of position and velocity vector to return
    vec Fx = zeros(nBodies);                // Gravitational force in x-direction to compute vxdot
    vec Fy = zeros(nBodies);                // Gravitational force in y-direction to compute vydot

    calculateForcesAndEnergy(posVel, Fx, Fy);

    // Finding the derivative of each velocity from the gravitational force
    for(int i = 0; i < nBodies; i++){
        dotPosVel[5*i+0] = posVel[5*i+2];               // setting xdot  = vx
        dotPosVel[5*i+1] = posVel[5*i+3];               // setting ydot  = vy
        dotPosVel[5*i+2] = Fx[i]/masses[i];             // setting vxdot = Fx/m
        dotPosVel[5*i+3] = Fy[i]/masses[i];             // setting xydot = Fy/m
    }   // Ending for-loop computing derivatives

    return dotPosVel;
}   // End of diffEq-function
*/
// FUNCTION TO USE IN CASE YOU HAVE TO COMMENT AWAY DIFFEQ

vec System::diffEq(vec, ofstream*){
    vec hei = vec(3);
    return hei;
}



void System::RK4(double h, ofstream* posvel, ofstream* enmom){
    // Function performing the Runge-Kutta 4 method

    // Declaring slopes to perform RK4
    vec k0, k1, k2, k3, k4;

    k0 = zeros(5*nBodies);
    k1 = h*diffEq(k0,    enmom);            //  IS THERE A GOOD WAY OF CREATING AND OPENING A FILE IN THE   //
    k2 = h*diffEq(k1/2., enmom);            //  CLASS ITSELF, SO THAT WE DON'T HAVE TO TAKE IT IN EVERY     //
    k3 = h*diffEq(k2/2., enmom);            //  SINGLE TIME?                                                //
    k4 = h*diffEq(k3,    enmom);

    // Updating position and velocity after RK4
    defineSystem = defineSystem + (k1 + 2.*k2 + 2.*k3 + k4)/6.;

    // Writing position and velocity to file
    *posvel << setiosflags(ios::showpoint | ios:: uppercase);
    *posvel << setw(20) << setprecision(15) << defineSystem.t();

}   // Ending RK4-function

/*
void calculateForcesAndEnergy(vec posVel, vec &Fx, vec &Fy){
    // Function calculating forces and energy (and angular momentum!) for the system

    double G = 4*M_PI*M_PI;                 // Defining the gravitational constant in appropriate units
    vec energyAngMom = zeros(4);            // Kinetic, potential and total energy. Angular momentum.
    Fx.zeros();
    Fy.zeros();

    for(int i = 0; i < nBodies; i++){
        for(int j=i+1; j < nBodies; j++){
            // Variables simplifying the calculations
            double dx     = posVel[5*i+0]-posVel[5*j+0];    // x-separation between celestial bodies
            double dy     = posVel[5*i+1]-posVel[5*j+1];    // y-separation between celestial bodies
            double dr2    = dx*dx + dy*dy;
            double dr     = sqrt(dr2);                      // separation between celestial bodies
            double factor = G*masses(j)*masses(i) / pow(dr,3);

            // Updating gravitational force and potential energy experienced by celestial object
            energyAngMom(1) += factor*dr;                   // Potential energy
            Fx[i] -= factor*dx;                             // Newton's law of gravity
            Fy[i] -= factor*dy;                             // Newton's law of gravity
            Fx[j] += factor*dx;                             // Using N3 to find the force on the opposite object
            Fy[j] += factor*dy;                             // Using N3 to find the force on the opposite object
        }   // Ending for-loop computing force and potential energy

        // Variables simplifying the calculations
        double pos = sqrt(pow(posVel[5*i+0],2) + pow(posVel[5*i+1],2)); // Radial vector/position in CoM frame
        double v2  = pow(posVel[5*i+2],2) + pow(posVel[5*i+3],2);       // Velocity squared of celestial object

        // Updating kinetic energy and angular momentum
        energyAngMom(3) += masses(i)*pos*sqrt(v2);                      // Angular momentum
        energyAngMom(0) += 0.5*masses(i)*v2;                            // Kinetic energy

    }   // Ending for-loop going over all celestial bodies

    energyAngMom(2) = sqrt(pow(energyAngMom(0),2) + pow(energyAngMom(1),2));   // Total energy

    // Writing energy and angular momentum to file
    *enmom << setiosflags(ios::showpoint | ios:: uppercase);
    *enmom << setw(20) << setprecision(15) << energyAngMom.t();

}   // Ending calculateForcesAndEnergy-function
*/








// ================================= MAIN FUNCTION ============================================== //
int mainInClasses(){

    double T = 50;                  // Total time of simulation (yrs)
    double h = 1/12;                // Step size
    int n = T/h;                    // Number of iterations
    int nBodies = 2;


    return 0;
}
