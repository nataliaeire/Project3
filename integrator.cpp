#include "integrator.h"
#include "cpelapsedtimer.h"
#include <stdlib.h>
using std::cout;
using std::endl;

Integrator::Integrator()
{ // Constructing integrator
    numThreads       = 1;
    counter          = 0; // Variable for accessing VerletInitialise
    adaptive_counter = 0;
    n                = 0;
} // End of constructor

Integrator::Integrator(int numThreads)
{ // Constructing integrator
    this->numThreads = numThreads;
    counter          = 0; // Variable for accessing VerletInitialise
    adaptive_counter = 0;
    n                = 0;
} // End of constructor


// ============================================= RK4 ==================================================== /
void Integrator::RK4(System &system, double dt)
{ // Function performing RK4 integration (by calling a new function to evolve systems inside)
    System k1 = system;
    k1.calculateForcesAndEnergy();
    //findDerivatives(k1); // f(t, y)

    System k2 = system;
    evolveSystem1InTimeUsingDerivativesFromSystem2(k2, k1, dt/2); // r = changed using v from k1
    k2.calculateForcesAndEnergy();
    //findDerivatives(k2); // f(t+dt/2, y + k1*dt/2)

    System k3 = system;
    evolveSystem1InTimeUsingDerivativesFromSystem2(k3, k2, dt/2); // r = changed using v from k2
    k3.calculateForcesAndEnergy();
    //findDerivatives(k3); // f(t+dt/2, y + k2*dt/2)

    System k4 = system;
    evolveSystem1InTimeUsingDerivativesFromSystem2(k4, k3, dt);   // r = changed using v from k3
    k4.calculateForcesAndEnergy();
    //findDerivatives(k4); // f(t+dt/2, y + k3*dt)

    for(int i=0; i<system.numberOfBodies(); i++){
        CelestialBody &body = system.bodies[i];

        double dt6 = dt/6;
        body.position.addAndMultiply(k1.bodies[i].velocity, dt6);
        body.position.addAndMultiply(k2.bodies[i].velocity, 2*dt6);
        body.position.addAndMultiply(k3.bodies[i].velocity, 2*dt6);
        body.position.addAndMultiply(k4.bodies[i].velocity, dt6);

        double dt6mass = dt6/body.mass;
        body.velocity.addAndMultiply(k1.bodies[i].force, dt6mass);
        body.velocity.addAndMultiply(k2.bodies[i].force, 2*dt6mass);
        body.velocity.addAndMultiply(k3.bodies[i].force, 2*dt6mass);
        body.velocity.addAndMultiply(k4.bodies[i].force, dt6mass);

    } // End for-loop

    system.calculateForcesAndEnergy();
} // End of RK4-function


void Integrator::RK4GR(System &system, double dt)
{ // Function performing RK4 integration for GR case (by calling a new function to evolve systems inside)
    System k1 = system;
    k1.calculateForcesUsingGR();
    //findDerivatives(k1); // f(t, y)

    System k2 = system;
    evolveSystem1InTimeUsingDerivativesFromSystem2(k2, k1, dt/2); // r = changed using v from k1
    k2.calculateForcesUsingGR();
    //findDerivatives(k2); // f(t+dt/2, y + k1*dt/2)

    System k3 = system;
    evolveSystem1InTimeUsingDerivativesFromSystem2(k3, k2, dt/2); // r = changed using v from k2
    k3.calculateForcesUsingGR();
    //findDerivatives(k3); // f(t+dt/2, y + k2*dt/2)

    System k4 = system;
    evolveSystem1InTimeUsingDerivativesFromSystem2(k4, k3, dt);   // r = changed using v from k3
    k4.calculateForcesUsingGR();
    //findDerivatives(k4); // f(t+dt/2, y + k3*dt)

    for(int i=0; i<system.numberOfBodies(); i++){
        CelestialBody &body = system.bodies[i];
        double dt6 = dt/6;
        body.position.addAndMultiply(k1.bodies[i].velocity, dt6);
        body.position.addAndMultiply(k2.bodies[i].velocity, 2*dt6);
        body.position.addAndMultiply(k3.bodies[i].velocity, 2*dt6);
        body.position.addAndMultiply(k4.bodies[i].velocity, dt6);

        double dt6mass = dt6/body.mass;
        body.velocity.addAndMultiply(k1.bodies[i].force, dt6mass);
        body.velocity.addAndMultiply(k2.bodies[i].force, 2*dt6mass);
        body.velocity.addAndMultiply(k3.bodies[i].force, 2*dt6mass);
        body.velocity.addAndMultiply(k4.bodies[i].force, dt6mass);
    } // End for-loop


    system.calculateForcesUsingGR();
} // End of RK4GR-function


void Integrator::evolveSystem1InTimeUsingDerivativesFromSystem2(System &system1, System &system2, double dt)
{ // Function evolving a system in time using values at previous time step
    for(int i=0; i<system1.numberOfBodies(); i++){
        CelestialBody &body1 = system1.bodies[i];
        CelestialBody &body2 = system2.bodies[i];
        body1.position.addAndMultiply(body2.velocity, dt);
        body1.velocity.addAndMultiply(body2.force, dt/body2.mass);
    } // End for-loop
} // End of evolveSystem1InTimeUsingDerivativesFromSystem2-function


// ============================================= VERLET ===================================================== //
void Integrator::VerletInitialise(System &system, double dt)
{ // Creates the system at time -dt. This system is used the first time the Verlet algorithm is run
    oldSystem = system;
    oldSystem.calculateForcesAndEnergy();
    for(int i=0; i < oldSystem.numberOfBodies(); i++){  // Loop changing all positions to what they were at time -dt
        CelestialBody &body1 = oldSystem.bodies[i];
        body1.position.addAndMultiply(body1.velocity, -dt);
        body1.position.addAndMultiply(body1.force, 0.5*dt*dt/body1.mass);
    } // End for-loop

    // Changing counter to avoid calling VerletInitialise next time Verlet is run
    counter = 1;
} // End of VerletInitialise-function


void Integrator::Verlet(System &system, double dt)
{ // Function performing Verlet integration
    // The first time Verlet() is run, VerletInitialise() is called
    if (counter == 0)   VerletInitialise(system, dt);

    System nextOldSystem = system;                      // Copying the system at time t to update old_system later on
    system.calculateForcesAndEnergy();                  // Calculates the forces on the bodies
    VerletEvolve(system, dt);                           // Evolving the system according to the Verlet algorithm
    oldSystem = nextOldSystem;                          // Updating oldSystem
} // End of Verlet-function


void Integrator::VerletEvolve(System &system, double dt)
{ // Function evolving the system for the Verlet integration
    for(int i=0; i < system.numberOfBodies(); i++){     // Looping over all bodies
        CelestialBody &body = system.bodies[i];         // the body at time t
        CelestialBody &bodyOld = oldSystem.bodies[i];   // the body at time t-dt

        vec3 tempPosDiff = body.position-bodyOld.position;
        body.position.add(tempPosDiff);
        body.position.addAndMultiply(body.force, dt*dt/body.mass);  // body.position is now at t+dt
        body.velocity = (body.position - bodyOld.position)/(2*dt);  // Calculating the velocity
    } // Ending for-loop
} // Ending VerletEvolve-function



// ======================================= VELOCITY VERLET ================================================== //
void Integrator::VelocityVerlet(System &system, double dt)
{ // Function performing Velocity Verlet integration

    // Initialisation
    if(counter == 0){
        system.calculateForcesAndEnergy();            // Calculates the forces on the bodies
        counter = 1;                                  // Updating the counter so that the if test is false
    } // End if-statement

    VelocityVerletEvolve(system, dt);                 // Evolving the system according to the Verlet algorithm

} // End of VelocityVerlet-function


void Integrator::VelocityVerletEvolve(System &system, double dt)
{ // Function evolving the system for the Velocity Verlet integration
//#pragma omp parallel for num_threads(numThreads)
    for(int i=0; i < system.numberOfBodies(); i++){     // Looping over all bodies
        CelestialBody &body = system.bodies[i];         // the body at time t

        // Calculating the velocity
        body.velocity.addAndMultiply(body.force, 0.5*dt/body.mass);
        body.position.addAndMultiply(body.velocity, dt);            // Calculating the position
    } // Ending for-loop

    system.calculateForcesAndEnergy();                  // Calculating the forces with new positions

    for(int i=0; i < system.numberOfBodies(); i++){     // Looping over all bodies
        CelestialBody &body = system.bodies[i];         // the body at time t+dt (velocity at time t)
        body.velocity.addAndMultiply(body.force, 0.5*dt/body.mass); // Calculating the velocity
    } // Ending for-loop
} // Ending VerletEvolve-function


// ====================================== ADAPTIVE VELOCITY VERLET =========================================== //
double Integrator::adaptiveDt()
{ // Simply returning the time step
    return adaptive_dt;
} // Ending adaptiveDt-function


void Integrator::initialiseAdaptiveVelocityVerlet(System &system)
{ // Function for initialising Velocity Verlet, as the forces need to be calculated
    for(int i = 0; i < int(system.numberOfBodies()); i++){
        CelestialBody &body = system.bodies[i];
        body.resetForce();
        system.actuallyCalculatingForces(body, n);
    } // End for-loop
} // End of initialiseAdaptiveVelocityVerlet-function


void Integrator::adaptiveVelocityVerlet(System &system)
{ // Function performing Velocity Verlet integration with adaptive time steps
    vec3 bodyacc;              // Variable for easily accessing the acceleration vector of a body
    double temp_acceleration;  // Variable for easily accessing the acceleration of a body
    double max_acc;            // Variable for storing the largest acceleration

    // Initialise for the very first time step (forces need to be calculated before starting
    if(adaptive_counter == 0)   initialiseAdaptiveVelocityVerlet(system);

    // Sort bodies and change time step every 1000 steps
    if(adaptive_counter % int(5) == 0){
        max_acc = 1e-7;                 // Sets an initial, very low acceleration
        CPElapsedTimer::sortBodies().start();
        system.sortBodiesIntoGroups();  // Updates bodies1, bodies2, etc.
        CPElapsedTimer::sortBodies().stop();

        // Finding the smallest time step:
        for(int j=0; j<int(system.bodies1.size()); j++){ // Finding the maximum acceleration in the system
            CelestialBody *body = system.bodies1[j];
            bodyacc = body->acceleration();
            temp_acceleration = bodyacc.length();
            if(temp_acceleration > max_acc) max_acc = temp_acceleration;
        } // End for-loop

        adaptive_dt = 1/max_acc;                    // Finds the time-step from the smallest acceleration in the system
        if(adaptive_dt < 1e-6) adaptive_dt = 1e-6;  // Overwrites the smallest time step if the time step is too low
    } // End if-statement

    system.resetEnergy();                           // Reset energy of system

    // The actual computations of Velocity Verlet
    for(n=0; n < 8; n++){
        halfKickAdaptively(system);
        moveBodies(system);
        fullKick(system);
    } // End foor-loop

    adaptive_counter++;

} // End adaptiveVelocityVerlet


void Integrator::halfKickAdaptively(System &system)
{ // Function calculating forces and energy (and angular momentum!) for the system

    double dt;
    if(n % 8 == 0){
        dt = 8*adaptive_dt;
        halfKick(system.bodies4, dt);
    } // Ending if-statement

    if(n % 4 == 0){
        dt = 4*adaptive_dt;
        halfKick(system.bodies3, dt);
    } // Ending if-statement

    if(n % 2 == 0){
        dt = 2*adaptive_dt;
        halfKick(system.bodies2, dt);
    } // Ending if-statement

    dt = adaptive_dt;
    halfKick(system.bodies1, dt);

} // Ending calculateForcesAndEnergy-function


void Integrator::calculateForcesForGroup(System &system, std::vector<CelestialBody*> &bodies)
{ // Function to calculate the forces between a specific group and all bodies
    int size = int(bodies.size());

    CPElapsedTimer::calculateForces().start();
#pragma omp parallel for shared(system, bodies, size) num_threads(numThreads)
    for(int i = 0; i < size; i++){
        CelestialBody *body = bodies[i];
        body->resetForce();
        system.actuallyCalculatingForces(*body, n);
    } // Ending for-loop
    CPElapsedTimer::calculateForces().stop();

    CPElapsedTimer::gatherEnergies().start();
    if(n == 7) system.gatherEnergyFromBodiesInGroup(bodies);
    CPElapsedTimer::gatherEnergies().stop();

} // End of calculateForcesForGroup-function


void Integrator::moveBodies(System &system)
{ // Moving bodies according to their velocity (2nd step in Velocity Verlet)
    for(int i=0; i<int(system.bodies.size()); i++){
        CelestialBody &body = system.bodies[i];
        body.position.addAndMultiply(body.velocity,adaptive_dt);
    } // Ending for-loop
} // Ending moveBodies-function


void Integrator::fullKick(System &system)
{ // Calculating forces and computing velocity at final time step (by calling halfKick)
    double dt;

    if((n+1) % 8 == 0){
        dt = 8*adaptive_dt;
        calculateForcesForGroup(system, system.bodies4);
        halfKick(system.bodies4, dt);
    } // Ending if-statement

    if((n+1) % 4 == 0){
        dt = 4*adaptive_dt;
        calculateForcesForGroup(system, system.bodies3);
        halfKick(system.bodies3, dt);
    } // Ending if-statement

    if((n+1) % 2 == 0){
        dt = 2*adaptive_dt;
        calculateForcesForGroup(system, system.bodies2);
        halfKick(system.bodies2, dt);
    } // Ending if-statement

    dt = adaptive_dt;
    calculateForcesForGroup(system, system.bodies1);
    halfKick(system.bodies1,dt);
} // Ending afterKick-function


void Integrator::halfKick(std::vector<CelestialBody*> &bodies, double dt)
{ // Calculating the velocity (step 1 and 3)
    for(int i=0; i<int(bodies.size()); i++){  // Looping over all bodies in group
        CelestialBody *body = bodies[i];
        body->velocity = body->velocity + body->acceleration()*dt*0.5;
    } // End for-loop
} // End halfKick-function

