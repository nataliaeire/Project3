#ifndef RK4_H
#define RK4_H
#include <system.h>
#include <celestialbody.h>

class Integrator
{
private:
    System old_system;
    int counter;        // Variable for accessing VerletInitialise
    void evolveSystem1InTimeUsingDerivativesFromSystem2(System &system1, System &system2, double dt);
    void VerletInitialise(System &system, double dt);
    void VerletEvolve(System &system1, double dt);


public:
    Integrator();
    void RK4(System &system, double dt);
    void Verlet(System &system, double dt);

};

#endif // RK4_H
