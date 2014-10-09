#include "system.h"
using namespace arma;

System::System() {
    angularMomentum = zeros(3);
    potentialEnergy = 0;
    kineticEnergy = 0;
}


void System::addBody(vec position, vec velocity, double mass){
    bodies.push_back( CelestialBody(position, velocity, mass) );
} // End addBody-function

void System::addBody(double x, double y, double z, double vx, double vy, double vz, double mass) {
    vec pos(3);
    vec vel(3);
    pos[0] = x; pos[1] = y; pos[2] = z;
    vel[0] = vx; vel[1] = vy; vel[2] = vz;
    addBody(pos, vel, mass);
} // End addBody-function

void System::calculateForcesAndEnergy(){
    // Function calculating forces and energy (and angular momentum!) for the system

    double G = 4*M_PI*M_PI;                 // Defining the gravitational constant in appropriate units
    potentialEnergy = 0;
    kineticEnergy = 0;
    angularMomentum.zeros();

    // Remember to reset forces before we calculate new ones
    for(int i=0; i<numberOfBodies(); i++) {
        CelestialBody &body = bodies[i];
        body.resetForce();
    }

    for(int i = 0; i < numberOfBodies(); i++){
        CelestialBody &body1 = bodies[i];

        for(int j=i+1; j < numberOfBodies(); j++){
            CelestialBody &body2 = bodies[j];
            vec deltaRVector = body1.position - body2.position;
            double dr = sqrt(dot(deltaRVector, deltaRVector));

            // Variables simplifying the calculations
            double factor = G*body1.mass*body2.mass / pow(dr,3);

            // Updating gravitational force and potential energy experienced by celestial object
            potentialEnergy += factor*dr;
            body1.force -= factor*deltaRVector;
            body2.force += factor*deltaRVector;
        }   // Ending for-loop computing force and potential energy

        // Variables simplifying the calculations
        kineticEnergy += 0.5*body1.mass*dot(body1.velocity, body1.velocity);
        vec momentum = body1.velocity*body1.mass;
        angularMomentum += cross(body1.position, momentum);
    }   // Ending for-loop going over all celestial bodies
}
// Ending calculateForcesAndEnergy-function

int System::numberOfBodies()
{
    return bodies.size();
}

double System::totalEnergy()
{
    return potentialEnergy + kineticEnergy;
}

