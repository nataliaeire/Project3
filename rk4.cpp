#include "rk4.h"

RK4::RK4()
{
}

void RK4::integrate(System &system, double dt)
{
    System k1 = system;
    findDerivatives(k1); // f(t, y)

    System k2 = system;
    evolveSystem1InTimeUsingDerivativesFromSystem2(k2, k1, dt/2); // r = changed using v from k1
    findDerivatives(k2); // f(t+dt/2, y + k1*dt/2)

    System k3 = system;
    evolveSystem1InTimeUsingDerivativesFromSystem2(k3, k2, dt/2); // r = changed using v from k2
    findDerivatives(k3); // f(t+dt/2, y + k2*dt/2)

    System k4 = system;
    evolveSystem1InTimeUsingDerivativesFromSystem2(k4, k3, dt); // r = changed using v from k3
    findDerivatives(k4); // f(t+dt/2, y + k3*dt)

    for(int i=0; i<system.numberOfBodies(); i++) {
        CelestialBody &body = system.bodies[i];
        body.position += 1.0/6*dt*(k1.bodies[i].positionDot + 2*k2.bodies[i].positionDot + 2*k3.bodies[i].positionDot + k4.bodies[i].positionDot);
        body.velocity += 1.0/6*dt*(k1.bodies[i].velocityDot + 2*k2.bodies[i].velocityDot + 2*k3.bodies[i].velocityDot + k4.bodies[i].velocityDot);
    }

    system.calculateForcesAndEnergy();
}

void RK4::findDerivatives(System &system)
{
    system.calculateForcesAndEnergy();

    for(int i=0; i<system.numberOfBodies(); i++) {
        CelestialBody &body = system.bodies[i];
        body.positionDot = body.velocity;
        body.velocityDot = body.force/body.mass;
    }
}

void RK4::evolveSystem1InTimeUsingDerivativesFromSystem2(System &system1, System &system2, double dt)
{
    for(int i=0; i<system1.numberOfBodies(); i++) {
        CelestialBody &body1 = system1.bodies[i];
        CelestialBody &body2 = system2.bodies[i];
        body1.position += body2.positionDot*dt;
        body1.velocity += body2.velocityDot*dt;
    }
}
