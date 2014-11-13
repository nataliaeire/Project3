#include "system.h"
#include <gaussiandeviate.h>
#include <cmath>
#include <armadillo>



// =================================== INITIALISING ================================== //
System::System()
{ // Initialising some qualities when creating a system, regardless of number of particles
    setG(false);
    potentialEnergy         = 0;
    kineticEnergy           = 0;
    boundKineticEnergy      = 0;
    boundPotentialEnergy    = 0;
    virialKineticEnergy     = 0;
    virialPotentialEnergy   = 0;
    angularMomentum.setToZero();
    totalMass               = 0;
    density                 = 0;
    smoothing               = false;
} // End constructor


void System::resetEnergy()
{ // Resetting all energies and angular momentum of system
    potentialEnergy         = 0;
    kineticEnergy           = 0;
    boundKineticEnergy      = 0;
    boundPotentialEnergy    = 0;
    virialKineticEnergy     = 0;
    virialPotentialEnergy   = 0;
    angularMomentum.setToZero();
} // End of resetEnergy-function


int System::numberOfBoundBodies()
{ // Finding number of bound celestial bodies in the system
    int boundBodies = 0;
    for(int i = 0; i < numberOfBodies(); i++){
        CelestialBody &body = bodies[i];
        if(body.gravitationallyBound == true)   boundBodies ++;
    } // Ending for-loop
    return boundBodies;
} // End of numberOfBoundBodies-function


void System::virialEnergy()
{ // Finding the virial energy of the system from the energy of bound bodies in the system
    virialKineticEnergy     = boundKineticEnergy/numberOfBoundBodies();
    virialPotentialEnergy   = boundPotentialEnergy/numberOfBoundBodies();
} // End of virialEnergy-function


void System::setG(bool cluster)
{ // Setting gravitational constant G based on the type of system
    if(cluster == 1){
        G = 3*M_PI / (32*density);              // t^2 = 3pi/32Grho, but t == 1
    }else{
        G = 4*M_PI*M_PI;
    } // End if-statement
} // End setG-function


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
    } // Ending for-loop
    CelestialBody &sun = bodies[0];
    sun.velocity = momentum/(-1*sun.mass);
} // End of conserveMomentum-function


void System::sortBodiesIntoGroups()
{ // Sorting bodies into groups based on their acceleration in order to perform adaptive steplength integration
    calculateForcesAndEnergy();
    arma::vec accelerations;
    accelerations.zeros(numberOfBodies());

    for(int i=0; i < numberOfBodies(); i++){
        CelestialBody &body = bodies[i];
        accelerations[i]    = log(body.acceleration().length());
    } // End for-loop

    // Sorting the accelerations
    double maxAcc   = arma::max(accelerations);
    double minAcc   = arma::min(accelerations);
    double deltaAcc = 0.25*(maxAcc - minAcc);

    // Initializing body groups
    bodies1.clear();
    bodies2.clear();
    bodies3.clear();
    bodies4.clear();

    // Looping over bodies to put them into separate groups
    for(int i = 0; i < numberOfBodies(); i++){
        CelestialBody &body  = bodies[i];
        if(log(body.acceleration().length()) < minAcc + deltaAcc){
            bodies4.push_back(&body);
        }else if(log(body.acceleration().length()) >= minAcc + deltaAcc && log(body.acceleration().length()) < minAcc + 2*deltaAcc){
            bodies3.push_back(&body);
        }else if(log(body.acceleration().length()) >= minAcc + 2*deltaAcc && log(body.acceleration().length()) < minAcc + 3*deltaAcc){
            bodies2.push_back(&body);
        }else{
            bodies1.push_back(&body);
        } // Ending if-statement
    } // Ending for-loop

} // End sortBodiesIntoGroups-function


// =================================== CREATING SYSTEM ================================== //
void System::addBody(vec3 position, vec3 velocity, double mass)
{ // Adding a celestial body, defined according to the CelestialBody-class to a system bodies
    bodies.push_back( CelestialBody(numberOfBodies(), position, velocity, mass) );
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


void System::addRandomSystem(int numberOfObjects, double sphereRadius)
{ // Adding system using random number generators
    double u, v, w, r, theta, phi, x, y, z, vx, vy, vz, massDeviation, mass;
    long int seed = 2;  // Seed to start random number generator
    totalMass = density = 0;

    // Loop to generate numberOfObjects celestial objects with random positions and mass
    for(int i = 0; i < numberOfObjects; i++){
        // Use that r^2 sin theta dr dtheta dphi = A du dv dw => r = R0(u)^(1/3), theta = arccos(1-2v), phi = 2pi w
        // x = r sin theta cos phi, y = r sin theta sin phi, z = r cos theta
        // Thus, generating random u, v & w will allow to get a uniform, spherical distribution

        // Generating random numbers (uniform distribution)
        u = ran2(&seed);
        v = ran2(&seed);
        w = ran2(&seed);

        // Calculating random spherical coordinates
        r       = sphereRadius*pow(u,1./3);
        theta   = acos(1.-2.*v);
        phi     = 2*M_PI*w;

        // Calculating random cartesian coordinates
        x = r*sin(theta)*cos(phi);
        y = r*sin(theta)*sin(phi);
        z = r*cos(theta);

        // Starting the system at rest
        vx = vy = vz = 0;

        // Generating random mass (normal distribution)
        massDeviation = gaussian_deviate(&seed);    // mean = 0, std = 1
        mass = 10. + massDeviation;                 // 10 solar masses + randomly drawn gaussian deviation
        totalMass += mass;                          // updating the total mass of the system

        // Adding body to system using previously created function
        addBody(x, y, z, vx, vy, vz, mass);
    } // Ending for-loop

    density = 3*totalMass / (4*M_PI*pow(sphereRadius,3));   // Density of the system
    setG(true);

} // End of addRandomSystem-function


// ===================================== PROPERTIES OF SYSTEM ========================================= //
int System::numberOfBodies()
{ // Simply returning the number of celestial bodies in the system
    return bodies.size();
} // End of numberOfBodies-function

double System::totalEnergy()
{ // Simply calculating the total energy of the system
    return potentialEnergy + kineticEnergy;
} // End of totalEnergy-system


// =================================== CALCULATING FORCES & ENERGY ===================================== //
void System::calculateForcesAndEnergy()
{ // Function calculating forces and energy (and angular momentum!) for the system
    // Initialising values
    resetEnergy();

    // Remembering to reset forces before we calculate new ones
    for(int i=0; i<numberOfBodies(); i++){
        CelestialBody &body = bodies[i];
        body.resetForce();
    } // Ending for-loop

    for(int i = 0; i < numberOfBodies(); i++){
        CelestialBody &body1 = bodies[i];

        for(int j=i+1; j < numberOfBodies(); j++){
            CelestialBody &body2 = bodies[j];

            // Variables simplifying the calculations
            vec3   deltaRVector = body1.position - body2.position;          // Spatial separation in all three directions
            double dr           = deltaRVector.length();                    // Separation radius/length/distance
            double factor       = G*body1.mass*body2.mass / pow(dr,3);      // Reoccuring factor

            // Updating the potential energy of the system and the gravitational force experienced by celestial object
            potentialEnergy -= factor*dr*dr;                                // Definition of the potential energy
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
    double c = 63239.7263;      // Defining the speed of light in units [AU/year]

    // Remembering to reset forces before we calculate new ones
    for(int i=0; i<numberOfBodies(); i++){
        CelestialBody &body = bodies[i];
        body.resetForce();
    } // Ending for-loop

    for(int i = 0; i < numberOfBodies(); i++){
        CelestialBody &body1 = bodies[i];

        for(int j=i+1; j < numberOfBodies(); j++){
            CelestialBody &body2 = bodies[j];

            // Variables simplifying the calculations
            vec3   deltaRVector     = body1.position - body2.position;                  // Spatial separation in 3D
            double dr               = deltaRVector.length();                            // Separation
            double NewtonianFactor  = G*body1.mass*body2.mass / pow(dr,3);              // Reoccuring factor
            vec3   angMomPerMass    = body2.position.cross(body2.velocity);             // L/m = r x v
            double GRFactor         = 1+(3*angMomPerMass.lengthSquared())/(dr*dr*c*c);  // Reoccuring factor
            double factor           = NewtonianFactor*GRFactor;

            // Updating gravitational force experienced by celestial object
            // Finding all components of the force
            body1.force.addAndMultiply(deltaRVector, -factor);                          // Law of gravity
            body2.force.addAndMultiply(deltaRVector, factor);                           // N3
        }   // Ending for-loop computing force and potential energy

    }   // Ending for-loop going over all celestial bodies

} // Ending calculateForcesAndEnergy-function


void System::actuallyCalculatingForces(CelestialBody &body, int n)
{ // Function finding the forces between
    // Initialising values
    body.resetEnergy();
    body.resetForce();

    for(int i = 0; i < numberOfBodies(); i++){
        CelestialBody &allBodies = bodies[i];

        // Only calculate forces between different celestial bodies
        if(body.index != allBodies.index){
            // Variables simplifying the calculations
            vec3   deltaRVector = allBodies.position - body.position;       // Spatial separation in all three directions
            double dr           = deltaRVector.length();                    // Separation radius/length/distance
            double factor;

            // Reoccuring factor
            if(smoothing == false){ // No smoothing
                factor = G*allBodies.mass*body.mass / pow(dr,3);
            }else { // Smoothing
                double epsilon = 0.1;
                factor = G*allBodies.mass*body.mass / ((dr*dr + epsilon*epsilon)*dr);
            } // Ending if-statement

            // Updating gravitational force and potential energy experienced by celestial object
            body.force.addAndMultiply(deltaRVector, factor);                // Gravitational law

            // Add contribution to potential energy if a time step has passed for all bodies
            if(n == 7)  body.PE -= factor*dr*dr;                            // Definition of the potential energy
        } // Ending if-statement

    }// Ending for-loop

    // Calculate kinetic energy and angular momentum if one large time step has passed
    if(n == 7){
        vec3 momentum   = body.velocity*body.mass;                          // p = m*v
        body.angMom     = body.position.cross(momentum);                    // L = r x p
        body.KE         = 0.5*body.mass*body.velocity.lengthSquared();      // KE = mv^2/2
    } // Ending if-statement

} // Ending actuallyCalculatingForces-function


void System::gatherEnergyFromBodiesInGroup(std::vector<CelestialBody*> &group)
{ // Calculating energies for bodies and for system
    for(int i = 0; i < int(group.size()); i++){
        CelestialBody *body = group[i];
        kineticEnergy      += body->KE;
        potentialEnergy    += 0.5*body->PE;
        angularMomentum.add(body->angMom);

        if(body->totalEnergyOfBody() <= 0.){
            body->gravitationallyBound  = true;
            boundKineticEnergy         += body->KE;
            boundPotentialEnergy       += 0.5*body->PE;
        }else{
            body->gravitationallyBound  = false;
        } // Ending if-statement
        virialEnergy();  // Calculating virial energy from bound energy
    } // Ending for-loop
} // End of calculateEnergy-function

