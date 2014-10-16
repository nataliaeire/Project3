#include "system.h"
#include <cmath>
#include <armadillo>

System::System()
{ // Initialising some qualities when creating a system, regardless of number of particles
    potentialEnergy = 0;
    kineticEnergy   = 0;
} // End constructor


void System::addBody(vec3 position, vec3 velocity, double mass)
{ // Adding a celestial body, defined according to the CelestialBody-class to a system bodies
    bodies.push_back( CelestialBody(position, velocity, mass) );
    //bodies.push_back adds a CelestialBody-object to the end of the object bodies //
} // End addBody-function


void System::addBody(double x, double y, double z, double vx, double vy, double vz, double mass)
{ // Alternative way of adding a body, which is more intuitive and easy from, for instance, a main function
  // using the former addBody-function to add a body to the system as before
    addBody(vec3(x,y,z), vec3(vx, vy, vz), mass);
} // End addBody-function


void System::addSystem(std::fstream &file)
{ // Alternative way of adding an entire system, reading in initials from file
    int valuesPerBody = 7;                          // For keeping track of which numbers belong to which object
    arma::vec values(valuesPerBody);                // Vector for storing values defining one body
    values.zeros(valuesPerBody);                    // Initialising vector

    double a;                                       // Input from file
    int i = 0;                                      // Counter variable

    while(file >> a){                               // Looping over all numbers in the file
        values(i) = a;                              // Storing details about a planet in a vector
        i++;
        if(i==valuesPerBody){                       // If-test for adding a body to the system
            addBody(values(0), values(1), values(2), values(3), values(4), values(5), values(6));
            i = 0;                                  // Re-initialising counter to start filling in a new vector
        } // End if-statement
    } // End while-loop
} // End addBody-function


int System::numberOfBodies()
{ // Simply returning the number of celestial bodies in the system
    return bodies.size();
} // End of numberOfBodies-function


void System::conserveMomentum()
{ // Function to find the momentum of the planets and changing the Sun's momentum to ensure
  // the total momentum of the system is conserved
    vec3 momentumTemp;
    momentumTemp.setToZero();
    momentum.setToZero();

    // Finding the total momentum of all bodies except the Sun
    for(int i=1; i<numberOfBodies(); i++){
        CelestialBody &body = bodies[i];
        momentumTemp = body.velocity*body.mass;
        momentum.add(momentumTemp);
    }

    CelestialBody &sun = bodies[0];
    sun.velocity = momentum/(-1*sun.mass);
}


void System::calculateForcesAndEnergy()
{ // Function calculating forces and energy (and angular momentum!) for the system

    // Initialising values
    double G = 4*M_PI*M_PI;     // Defining the gravitational constant in appropriate units
    potentialEnergy = 0;
    kineticEnergy   = 0;
    angularMomentum.setToZero();

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
            vec3   deltaRVector = body1.position - body2.position;          // Spatial separation in all three directions
            double dr           = deltaRVector.length();                    // Separation radius/length/distance
            double factor       = G*body1.mass*body2.mass / pow(dr,3);      // Reoccuring factor

            // Updating gravitational force and potential energy experienced by celestial object
            potentialEnergy += factor*dr;                                   // Definition of the potential energy
            // Finding all components of the force
            body1.force.addAndMultiply(deltaRVector, -factor);              // Law of gravity
            body2.force.addAndMultiply(deltaRVector, factor);               // N3
        }   // Ending for-loop computing force and potential energy

        // Variables simplifying the calculations
        vec3 momentum   = body1.velocity*body1.mass;                        // p = m*v
        angularMomentum.add(body1.position.cross(momentum));                // L = r x p, updated for each body
        kineticEnergy  += 0.5*body1.mass*body1.velocity.lengthSquared();    // k = mv^2/2, updated for each body
    }   // Ending for-loop going over all celestial bodies

} // Ending calculateForcesAndEnergy-function


void System::calculateForcesUsingGR()
{ // Function calculating forces for the system taking into account GR
    // Initialising values
    double G = 4*M_PI*M_PI;     // Defining the gravitational constant in appropriate units
    double c = 63239.7263;      // Defining the speed of light in units [AU/year]

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
            vec3   deltaRVector     = body1.position - body2.position;          // Spatial separation in all three directions
            double dr               = deltaRVector.length();                    // Separation radius/length/distance
            double NewtonianFactor  = G*body1.mass*body2.mass / pow(dr,3);      // Reoccuring factor
            vec3   momentum         = body1.velocity*body1.mass;                // p = m*v
            vec3   angMom           = body1.position.cross(momentum);           // L = r x p, updated for each body
            double GRFactor         = 1+(3*angMom.lengthSquared())/(dr*dr*c*c); // Reoccuring factor
            double factor           = NewtonianFactor*GRFactor;

            // Updating gravitational force experienced by celestial object
            // Finding all components of the force
            body1.force.addAndMultiply(deltaRVector, -factor);                  // Law of gravity
            body2.force.addAndMultiply(deltaRVector, factor);                   // N3
        }   // Ending for-loop computing force and potential energy

    }   // Ending for-loop going over all celestial bodies

} // Ending calculateForcesAndEnergy-function


double System::totalEnergy()
{ // Simply calculating the total energy of the system
    return potentialEnergy + kineticEnergy;
} // End of totalEnergy-system

