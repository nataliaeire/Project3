#ifndef RK4_H
#define RK4_H
#include <system.h>
#include <celestialbody.h>

class Integrator
{
private:
    System  oldSystem;
    int     counter;        // Variable for accessing VerletInitialise
    void    evolveSystem1InTimeUsingDerivativesFromSystem2(System &system1, System &system2, double dt);
    void    VerletInitialise(System &system, double dt);
    void    VerletEvolve(System &system1, double dt);
    void    VelocityVerletEvolve(System &system1, double dt);

public:
    Integrator();
    void    RK4(System &system, double dt);
    void    RK4GR(System &system, double dt);
    void    Verlet(System &system, double dt);
    void    VelocityVerlet(System &system, double dt);

};

#endif // RK4_H
