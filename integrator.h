#ifndef RK4_H
#define RK4_H
#include <system.h>
#include <celestialbody.h>

class Integrator
{
private:
    System  oldSystem;
    int     counter;        // Variable for accessing VerletInitialise
    int     n;              // Counter which choose which groups to calculate forces on
    double  adaptive_dt;    // Variable used by adaptiveVelocityVerlet
    void    evolveSystem1InTimeUsingDerivativesFromSystem2(System &system1, System &system2, double dt);
    void    VerletInitialise(System &system, double dt);
    void    VerletEvolve(System &system1, double dt);
    void    VelocityVerletEvolve(System &system1, double dt);
    void    evolveRightBodies(System &system, CelestialBody &body);
    void    adaptiveVelocityVerletEvolve(System &system);
    void    moveBodiesLinearly(CelestialBody &body);
    void    moveBodies(System &system);


public:
    Integrator();
    void    RK4(System &system, double dt);
    void    RK4GR(System &system, double dt);
    void    Verlet(System &system, double dt);
    void    VelocityVerlet(System &system, double dt);
    void    adaptiveVelocityVerlet(System &system, int i);
    double  adaptiveDt();

};

#endif // RK4_H
