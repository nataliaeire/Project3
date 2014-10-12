#include "integrator.h"
#include <stdlib.h>
Integrator::Integrator()
{
    counter = 0; // Variable for accessing VerletInitialise
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

void Integrator::evolveSystem1InTimeUsingDerivativesFromSystem2(System &system1, System &system2, double dt)
{
    for(int i=0; i<system1.numberOfBodies(); i++){
        CelestialBody &body1 = system1.bodies[i];
        CelestialBody &body2 = system2.bodies[i];
        body1.position.addAndMultiply(body2.velocity, dt);
        body1.velocity.addAndMultiply(body2.force, dt/body2.mass);
    }
}

void Integrator::VerletInitialise(System &system, double dt)// Creates the system at time -dt. This system
{                                                           // is used the first time the Verlet algorithm is run
    old_system = system;
    for(int i=0; i < old_system.numberOfBodies(); i++)      // Loop changing all positions to what they were
    {                                                       // at time -dt.
        CelestialBody &body1 = old_system.bodies[i];
        body1.position.addAndMultiply(body1.velocity, -dt);
    }
    counter = 1;                                            // The counter is changed to avoid calling
}                                                           // VerletInitialise the next time Verlet is run


void Integrator::Verlet(System &system, double dt)
{
    if (counter == 0)                                      // The first time Verlet() is run,
    {                                                      // VerletInitialise() is called
        VerletInitialise(system, dt);
    }
    System next_old_system = system;                       // Copying the system at time t to update old_system later on
    system.calculateForcesAndEnergy();                     // Calculates the forces on the bodies
    VerletEvolve(system, dt);                              // Evolving the system according to the Verlet algorithm
    old_system = next_old_system;                          // Updating old_system
}

void Integrator::VerletEvolve(System &system, double dt)
{
    for(int i=0; i < system.numberOfBodies(); i++)         // Looping over all bodies
    {
        CelestialBody &body = system.bodies[i];            // the body at time t
        CelestialBody &body_old = old_system.bodies[i];    // the body at time t-dt

        body.position.add(body.position-body_old.position);
        body.position.addAndMultiply(body.force, dt*dt/body.mass); // body.position is now at t+dt
        body.velocity.addAndMultiply(body.position - body_old.position, 1./(2*dt)); // Calculating the velocity

    }
}
