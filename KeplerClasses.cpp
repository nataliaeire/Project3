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
}

void System::addBody(vec position, vec velocity, double mass){
    CelestialBody body(position, velocity, mass);
    bodies[nBodies] = body;
    this->nBodies ++;
}

void System::definingBodies(){
    // Creating system described by a vector defineSystem
    // with celestial bodies defined according to the CelestialBody class

    // The vector defineSystem is created to contain all position, velocity and mass information about the system
    defineSystem = zeros(5*nBodies);

    for(int i = 0; i < nBodies; i++){             // Defining the system by position, velocity and mass
        defineSystem[5*i]   = bodies[i].position[0];
        defineSystem[5*i+1] = bodies[i].position[1];
        defineSystem[5*i+2] = bodies[i].velocity[0];
        defineSystem[5*i+3] = bodies[i].velocity[1];
        defineSystem[5*i+4] = bodies[i].mass;
    }   // End for-loop defining the system
}

