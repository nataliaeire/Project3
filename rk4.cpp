#include "rk4.h"

RK4::RK4()
{
}

void RK4::integrate(System &system, double dt)
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
    evolveSystem1InTimeUsingDerivativesFromSystem2(k4, k3, dt); // r = changed using v from k3
    k4.calculateForcesAndEnergy();
    //findDerivatives(k4); // f(t+dt/2, y + k3*dt)

    for(int i=0; i<system.numberOfBodies(); i++) {
        CelestialBody &body = system.bodies[i];
        body.position.addAndMultiply(k1.bodies[i].velocity + k2.bodies[i].velocity + k2.bodies[i].velocity + k3.bodies[i].velocity + k3.bodies[i].velocity + k4.bodies[i].velocity, dt/6);
        body.velocity.addAndMultiply(k1.bodies[i].force + k2.bodies[i].force + k2.bodies[i].force + k3.bodies[i].force + k3.bodies[i].force + k4.bodies[i].force, 1.0/6*dt/body.mass);
    }

    system.calculateForcesAndEnergy();
}

void RK4::evolveSystem1InTimeUsingDerivativesFromSystem2(System &system1, System &system2, double dt)
{
    for(int i=0; i<system1.numberOfBodies(); i++) {
        CelestialBody &body1 = system1.bodies[i];
        CelestialBody &body2 = system2.bodies[i];
        body1.position.addAndMultiply(body2.velocity, dt);
        body1.velocity.addAndMultiply(body2.force, dt/body2.mass);
    }
}
