#include "integrator.h"
#include <stdlib.h>
Integrator::Integrator()
{
    counter = 0; // Variable for accessing VerletInitialise
    n = 1;
}

void Integrator::RK4(System &system, double dt)
{
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
        body.position.addAndMultiply(k1.bodies[i].velocity + k2.bodies[i].velocity + k2.bodies[i].velocity + k3.bodies[i].velocity + k3.bodies[i].velocity + k4.bodies[i].velocity, dt/6);
        body.velocity.addAndMultiply(k1.bodies[i].force + k2.bodies[i].force + k2.bodies[i].force + k3.bodies[i].force + k3.bodies[i].force + k4.bodies[i].force, 1.0/6*dt/body.mass);
    }

    system.calculateForcesAndEnergy();
}


void Integrator::RK4GR(System &system, double dt)
{
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
        body.position.addAndMultiply(k1.bodies[i].velocity + k2.bodies[i].velocity + k2.bodies[i].velocity + k3.bodies[i].velocity + k3.bodies[i].velocity + k4.bodies[i].velocity, dt/6);
        body.velocity.addAndMultiply(k1.bodies[i].force + k2.bodies[i].force + k2.bodies[i].force + k3.bodies[i].force + k3.bodies[i].force + k4.bodies[i].force, 1.0/6*dt/body.mass);
    }

    system.calculateForcesUsingGR();
}


void Integrator::evolveSystem1InTimeUsingDerivativesFromSystem2(System &system1, System &system2, double dt)
{
    for(int i=0; i<system1.numberOfBodies(); i++){
        CelestialBody &body1 = system1.bodies[i];
        CelestialBody &body2 = system2.bodies[i];
        body1.position.addAndMultiply(body2.velocity, dt);
        body1.velocity.addAndMultiply(body2.force, dt/body2.mass);
    }
}

void Integrator::VerletInitialise(System &system, double dt)
{ // Creates the system at time -dt. This system is used the first time the Verlet algorithm is run
    oldSystem = system;
    for(int i=0; i < oldSystem.numberOfBodies(); i++){  // Loop changing all positions to what they were at time -dt
        CelestialBody &body1 = oldSystem.bodies[i];
        body1.position.addAndMultiply(body1.velocity, -dt);
    }
    // Changing counter to avoid calling VerletInitialise next time Verlet is run
    counter = 1;
} // End VerletInitialise-function

void Integrator::Verlet(System &system, double dt)
{
    // The first time Verlet() is run, VerletInitialise() is called
    if (counter == 0)   VerletInitialise(system, dt);

    System nextOldSystem = system;                      // Copying the system at time t to update old_system later on
    system.calculateForcesAndEnergy();                  // Calculates the forces on the bodies
    VerletEvolve(system, dt);                           // Evolving the system according to the Verlet algorithm
    oldSystem = nextOldSystem;                          // Updating oldSystem
}

void Integrator::VerletEvolve(System &system, double dt)
{
    for(int i=0; i < system.numberOfBodies(); i++){     // Looping over all bodies
        CelestialBody &body = system.bodies[i];         // the body at time t
        CelestialBody &bodyOld = oldSystem.bodies[i];   // the body at time t-dt

        body.position.add(body.position-bodyOld.position);
        body.position.addAndMultiply(body.force, dt*dt/body.mass);     // body.position is now at t+dt
        body.velocity = (body.position - bodyOld.position)*1./(2*dt);  // Calculating the velocity
    } // Ending for-loop
} // Ending VerletEvolve-function

void Integrator::VelocityVerlet(System &system, double dt)
{
    system.calculateForcesAndEnergy();                  // Calculates the forces on the bodies
    VelocityVerletEvolve(system, dt);                   // Evolving the system according to the Verlet algorithm
}

void Integrator::VelocityVerletEvolve(System &system, double dt)
{
    for(int i=0; i < system.numberOfBodies(); i++){     // Looping over all bodies
        CelestialBody &body = system.bodies[i];         // the body at time t
        vec3 velocity_dt_2;

        // Calculating the velocity
        velocity_dt_2 = body.velocity + body.force/body.mass*dt/2.;
        body.position = body.position + velocity_dt_2*dt;              // Calculating the position
        system.calculateForcesAndEnergy();
        body.velocity = velocity_dt_2 + body.force/body.mass*dt/2;  // Calculating the velocity
    } // Ending for-loop
} // Ending VerletEvolve-function


void Integrator::adaptiveVelocityVerlet(System &system, int i){
    vec3 bodyacc;
    double temp_acceleration;
    double max_acc;

    if(i % int(1e3) == 0){
        max_acc = 1e-7;
        system.sortBodiesIntoGroups();  // Updates bodies1, bodies2, etc.
        for(int j=0; j<system.bodies1.size();j++){
            CelestialBody *body = system.bodies1[j];
            bodyacc = body->acceleration;
            temp_acceleration = bodyacc.length();
            if(temp_acceleration > max_acc) max_acc = temp_acceleration;
        } // End for-loop

        adaptive_dt = 5/max_acc;
        if(adaptive_dt < 1e-4)     adaptive_dt = 1e-4;

    } // End if-statement

    system.calculateForcesAdaptively(n, smoothing);                 // Calculates the forces on the bodies
    moveBodies(system);
    adaptiveVelocityVerletEvolve(system);            // Evolving the system according to the Verlet algorithm
    if(n==8)  n=0;
    n += 1;
}


double Integrator::adaptiveDt()
{
    return adaptive_dt;
}


void Integrator::moveBodies(System &system)
{ // Function determining which bodies to move, but not calculate forces on

    if(n != 8){
        for(int i = 0; i < int(system.bodies4.size()); i++){
        CelestialBody *body4 = system.bodies4[i];
        moveBodiesLinearly(*body4);
        } // Ending for-loop
    } // Ending if-statement
    if(n != 8 && n != 4){
        for(int i = 0; i < int(system.bodies3.size()); i++){
        CelestialBody *body3 = system.bodies3[i];
        moveBodiesLinearly(*body3);
        } // Ending for-loop
    } // Ending if-statement
    if(n != 8 && n != 6 && n != 4 && n != 2){
        for(int i = 0; i < int(system.bodies2.size()); i++){
        CelestialBody *body2 = system.bodies2[i];
        moveBodiesLinearly(*body2);
        } // Ending for-loop
      } // Ending if-statement
    } // Ending calculateForcesAndEnergy-function


void Integrator::moveBodiesLinearly(CelestialBody &body)
{ // Moving bodies according to their velocity
    body.position.add(body.velocity*adaptive_dt);
} // Ending moveBodiesLinearly


void Integrator::adaptiveVelocityVerletEvolve(System &system){
    if(n % 8 == 0){
        for(int i = 0; i < int(system.bodies4.size()); i++){
        CelestialBody *body4 = system.bodies4[i];
        evolveRightBodies(system, *body4);
        } // Ending for-loop
    } // Ending if-statement
    if(n % 4 == 0){
        for(int i = 0; i < int(system.bodies3.size()); i++){
        CelestialBody *body3 = system.bodies3[i];
        evolveRightBodies(system, *body3);
        } // Ending for-loop
    } // Ending if-statement
    if(n % 2 == 0){
        for(int i = 0; i < int(system.bodies2.size()); i++){
        CelestialBody *body2 = system.bodies2[i];
        evolveRightBodies(system, *body2);
        } // Ending for-loop
    } // Ending if-statement
    for(int i = 0; i < int(system.bodies1.size()); i++){
        CelestialBody *body1 = system.bodies1[i];
        evolveRightBodies(system, *body1);
    } // Ending for-loop
} // Ending adaptiveVelocityVerletEvolve-function


void Integrator::evolveRightBodies(System &system, CelestialBody &body)
{ // Moving bodies according to their forces
    vec3 velocity_dt_2;

    // Calculating the velocity
    velocity_dt_2 = body.velocity + body.force/body.mass*adaptive_dt/2.;
    body.position = body.position + velocity_dt_2*adaptive_dt;              // Calculating the position
    system.calculateForcesAdaptively(n);
    body.velocity = velocity_dt_2 + body.force/body.mass*adaptive_dt/2;  // Calculating the velocity
 } // Ending VerletEvolve-function
