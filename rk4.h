#ifndef RK4_H
#define RK4_H
#include <system.h>
#include <celestialbody.h>

class RK4
{
private:
    void findDerivatives(System &system);
    void evolveSystem1InTimeUsingDerivativesFromSystem2(System &system1, System &system2, double dt);

public:
    RK4();
    void integrate(System &system, double dt);

};

#endif // RK4_H
