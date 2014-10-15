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


void System::addBody(std::fstream &file)
{ // Alternative way of adding a body, reading in initials from file
    // THIS FUNCTION NEEDS TO ACTUALLY BE WRITTEN //

    std::cout << "In addBody" << std::endl;       // Going swimmingly so far :)

    int valuesPerBody = 7;                            // For keeping track of which numbers belong to which object
    arma::vec values = arma::zeros(valuesPerBody);              // Vector for storing values belonging to one body

    std::cout << "In addBody" << std::endl;       // We made it here as well :)


    for (int j=0; j< 2*valuesPerBody-1; j++)       // This needs to be modified
    {
        std::cout << "In while-loop" << std::endl; // ...because this won't print
        for(int i=0; i<valuesPerBody; i++)         // Idea: Looping over values belongig to one body
        {                                          // so we won't get a lot of indices to keep track of
            file >> values(i);
            if(i==valuesPerBody-1)                 // If-test for printing one body
            {
                // Adding said body:
                addBody(values(0), values(1), values(2), values(3), values(4), values(5), values(6));
            }

        }
    }
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


double System::totalEnergy()
{ // Simply calculating the total energy of the system
    return potentialEnergy + kineticEnergy;
} // End of totalEnergy-system

