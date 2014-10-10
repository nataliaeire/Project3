#ifndef RK4_H
#define RK4_H
#include <system.h>
#include <celestialbody.h>

class Integrator
{
private:
    void evolveSystem1InTimeUsingDerivativesFromSystem2(System &system1, System &system2, double dt);

public:
    Integrator();
    void RK4(System &system, double dt);
};

#endif // RK4_H
