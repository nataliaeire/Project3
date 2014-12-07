#ifndef RK4_H
#define RK4_H
#include <system.h>
#include <celestialbody.h>

class Integrator
{
private:
    System  oldSystem;        // The system at the previous time step, used by Verlet
    int     counter;          // Variable for accessing VerletInitialise
    int     adaptive_counter; // Counter to reorder the groups of bodies every now and then
    int     n;                // Counter to choose which groups to calculate forces on
    double  adaptive_dt;      // Variable used by adaptiveVelocityVerlet

    // RK4
    void    evolveSystem1InTimeUsingDerivativesFromSystem2(System &system1, System &system2, double dt);

    // Verlet
    void    VerletInitialise(System &system, double dt);
    void    VerletEvolve(System &system1, double dt);

    // Velocity Verlet
    void    VelocityVerletEvolve(System &system1, double dt);

    // Adaptive Velocity Verlet
    void    initialiseAdaptiveVelocityVerlet(System &system);
    void    moveBodies(System &system);
    void    halfKickAdaptively(System &system);
    void    fullKick(System &system);
    void    halfKick(std::vector<CelestialBody*> &bodies, double dt);
    void    calculateForcesForGroup(System &system, std::vector<CelestialBody*> &bodies);

public:
    Integrator();
    Integrator(int numThreads);  // Constructor
    int     numThreads;

    // Integrators called in main
    void    RK4(System &system, double dt);             // Runge-Kutta 4
    void    RK4GR(System &system, double dt);           // Runge-Kutta 4 with a correction from general relativity
    void    Verlet(System &system, double dt);          // Verlet
    void    VelocityVerlet(System &system, double dt);  // Velocity Verlet
    void    adaptiveVelocityVerlet(System &system);     // Velocity Verlet with an adaptive time step

    double  adaptiveDt();                               // Function returning the value of adaptive_dt (a private variable)

};  // End of Integrator class declaration

#endif // RK4_H
