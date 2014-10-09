#include "system.h"
using namespace arma;

System::System()
{ // Initialising some qualities when creating a system, regardless of number of particles
    angularMomentum = zeros(3);
    potentialEnergy = 0;
    kineticEnergy   = 0;
} // End constructor


void System::addBody(vec position, vec velocity, double mass)
{ // Adding a celestial body, defined according to the CelestialBody-class to a system bodies
    bodies.push_back( CelestialBody(position, velocity, mass) );
    //bodies.push_back adds a CelestialBody-object to the end of the object bodies //
} // End addBody-function


void System::addBody(double x, double y, double z, double vx, double vy, double vz, double mass)
{ // Alternative way of adding a body, which is more intuitive and easy from, for instance, a main function

    // Creating position and velocity vectors from the input
    vec pos(3);
    vec vel(3);
    pos[0] = x;     pos[1] = y;     pos[2] = z;
    vel[0] = vx;    vel[1] = vy;    vel[2] = vz;

    // Using the former addBody-function to add a body to the system as before
    addBody(pos, vel, mass);

} // End addBody-function


int System::numberOfBodies()
{ // Simply returning the number of celestial bodies in the system
    return bodies.size();
} // End of numberOfBodies-function


void System::calculateForcesAndEnergy()
{ // Function calculating forces and energy (and angular momentum!) for the system

    // Initialising values
    double G = 4*M_PI*M_PI;     // Defining the gravitational constant in appropriate units
    potentialEnergy = 0;
    kineticEnergy   = 0;
    angularMomentum.zeros();

    // Remembering to reset forces before we calculate new ones
    for(int i=0; i<numberOfBodies(); i++){
        CelestialBody &body = bodies[i];
        body.resetForce();
    }

    for(int i = 0; i < numberOfBodies(); i++){
        CelestialBody &body1 = bodies[i];

        for(int j=i+1; j < numberOfBodies(); j++){
            CelestialBody &body2 = bodies[j];

            // Variables simplifying the calculations
            vec    deltaRVector = body1.position - body2.position;          // Spatial separation in all three directions
            double dr           = sqrt(dot(deltaRVector, deltaRVector));    // Separation radius/length/distance
            double factor       = G*body1.mass*body2.mass / pow(dr,3);      // Reoccuring factor

            // Updating gravitational force and potential energy experienced by celestial object
            potentialEnergy += factor*dr;           // Definition of the potential energy
            body1.force     -= factor*deltaRVector; // Finding all components of the force using the law of gravity
            body2.force     += factor*deltaRVector; // Finding all components of the force on opposite object using N3
        }   // Ending for-loop computing force and potential energy

        // Variables simplifying the calculations
        vec momentum     = body1.mass*body1.velocity;                           // p = m*v
        angularMomentum += cross(body1.position, momentum);                     // L = r x p, updated for each body
        kineticEnergy   += 0.5*body1.mass*dot(body1.velocity, body1.velocity);  // k = mv^2/2, updated for each body
    }   // Ending for-loop going over all celestial bodies

} // Ending calculateForcesAndEnergy-function


double System::totalEnergy()
{ // Simply calculating the total energy of the system
    return potentialEnergy + kineticEnergy;
} // End of totalEnergy-system

