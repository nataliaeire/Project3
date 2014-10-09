#include "KeplerClasses.h"


// ============================================== CLASS: CELESTIALBODY ============================================== //

// Initialising functions
CelestialBody::CelestialBody(){
}

CelestialBody::CelestialBody(double x, double y, double z, double vx, double vy, double vz, double mass){
    // Function defines a celestial body by position, velocity and mass from the input parameters

    this->x    = x;
    this->y    = y;
    this->z    = z;
    this->vx   = vx;
    this->vy   = vy;
    this->vz   = vz;
    this->mass = mass;
}


// ============================================== CLASS: SYSTEM ===================================================== //

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
