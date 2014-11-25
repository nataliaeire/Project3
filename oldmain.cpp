
#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace arma;

vec diffEq(vec positionAndVelocity);
vec verletInitialise(vec positionAndVelocity, double dt);

void rk4(double h, vec &positionAndVelocity, ofstream *file);
vec verlet(vec &positionAndVelocity, vec oldPositionAndVelocity, double dt, ofstream *file);
void energyAndAngMom(vec &positionAndVelocity, ofstream* file);

int oldmain()
{
    int n = 1e6;                   // Number of iterations
    double T = 1;                   // Total time of simulation (yrs)
    double h = T/n;                 // Step size

    // Create vector containting all positions and velocities
    vec positionAndVelocity(4);
    positionAndVelocity(0) = 1;         // Initial x-position
    positionAndVelocity(1) = 0;         // Initial y-position
    positionAndVelocity(2) = 0;         // Initial vx-velocity
    positionAndVelocity(3) = 2*M_PI;    // Initial vy-velocity (2*1 AU*pi / 1 yr)

    vec positionAndVelocityVerlet(4);
    positionAndVelocityVerlet(0) = 1;       // Initial x-position
    positionAndVelocityVerlet(1) = 0;       // Initial y-position
    positionAndVelocityVerlet(2) = 0;       // Initial vx-velocity
    positionAndVelocityVerlet(3) = 2*M_PI;  // Initial vy-velocity (2*1 AU*pi / 1 yr)


    // Creating file to write to
    ofstream posVelOut;
    ofstream energyOut;
    ofstream posVelOutVerlet;
    ofstream energyOutVerlet;
    posVelOut.open("posvel_RK4_sunfixed_dte-6.txt");
    energyOut.open("energyangmom_RK4_sunfixeddte-6.txt");
    posVelOutVerlet.open("posvel_Verlet_sunfixed_dte-6.txt");
    energyOutVerlet.open("energyangmom_Verlet_sunfixeddte-6.txt");

    // Writing initial positions and velocities to file
    posVelOut << setiosflags(ios::showpoint | ios:: uppercase);
    posVelOut << setw(20) << setprecision(15) << positionAndVelocity.t();
    posVelOutVerlet << setw(20) << setprecision(15) << positionAndVelocityVerlet.t();


    // Getting the system at the previous time step for Verlet
    vec oldPositionAndVelocity = verletInitialise(positionAndVelocityVerlet, h);
    vec oldPositionAndVel;

    for(int i=0; i<n; i++){
        energyAndAngMom(positionAndVelocity, &energyOut);              // Running RK4
        rk4(h, positionAndVelocity, &posVelOut);

        energyAndAngMom(positionAndVelocityVerlet, &energyOutVerlet);  // Running Verlet
        oldPositionAndVel = verlet(positionAndVelocityVerlet, oldPositionAndVelocity, h, &posVelOutVerlet);
        oldPositionAndVelocity = oldPositionAndVel;

    }

    posVelOut.close();
    energyOut.close();

    return 0;
}   // End of main function

// Function performing the Runge-Kutta 4 method
void rk4(double h, vec &positionAndVelocity, ofstream* file)
{
    // Declaring slopes to perform RK4
    vec k1, k2, k3, k4;

    k1 = h*diffEq(positionAndVelocity);
    k2 = h*diffEq(positionAndVelocity+k1/2.);
    k3 = h*diffEq(positionAndVelocity+k2/2.);
    k4 = h*diffEq(positionAndVelocity+k3);

    positionAndVelocity = positionAndVelocity + (k1 + 2.*k2 + 2.*k3 + k4)/6.;

    // Writing position and velocity to file
    *file << setiosflags(ios::showpoint | ios:: uppercase);
    *file << setw(20) << setprecision(15) << positionAndVelocity.t();
}


vec verletInitialise(vec positionAndVelocity, double dt)
{   // Using relation between position, velocity and acceleration when the acceleration is constant
    // to approximate the system at a time t-h
    vec oldPositionAndVelocity = positionAndVelocity;

    // Using the velocities in x- and y-direction
    oldPositionAndVelocity(0) -= positionAndVelocity(2)*dt;
    oldPositionAndVelocity(1) -= positionAndVelocity(3)*dt;

    vec velacc = diffEq(positionAndVelocity);  // Getting the acceleration

    // Using the accelerations in x- and y-direction
    oldPositionAndVelocity(0) += 0.5*dt*dt*velacc(2);
    oldPositionAndVelocity(1) += 0.5*dt*dt*velacc(3);


    return oldPositionAndVelocity;
}

vec verlet(vec &positionAndVelocity, vec oldPositionAndVelocity, double dt, ofstream *file)
{
    vec next_old_posvel = positionAndVelocity;

    vec velacc = diffEq(positionAndVelocity);

    positionAndVelocity(0) += positionAndVelocity(0)- oldPositionAndVelocity(0) + velacc(2)*dt*dt;
    positionAndVelocity(1) += positionAndVelocity(1)- oldPositionAndVelocity(1) + velacc(3)*dt*dt;
    positionAndVelocity(2) = (positionAndVelocity(0) - oldPositionAndVelocity(0))/(2.0*dt);
    positionAndVelocity(3) = (positionAndVelocity(1) - oldPositionAndVelocity(1))/(2.0*dt);
    oldPositionAndVelocity = next_old_posvel;

    // Writing position and velocity to file
    *file << setiosflags(ios::showpoint | ios:: uppercase);
    *file << setw(20) << setprecision(15) << positionAndVelocity.t();

   return oldPositionAndVelocity;
}

// Function setting the differential equations
vec diffEq(vec positionAndVelocity)
{
    vec posVel(4);      // Creating temporary position and velocity vector
    double r = sqrt( pow(positionAndVelocity[0],2)+pow(positionAndVelocity[1],2) );   // Calculating radius

    posVel[0] = positionAndVelocity[2];                          // setting xdot = vx
    posVel[1] = positionAndVelocity[3];                          // setting ydot = vy
    posVel[2] = -4.*M_PI*M_PI*positionAndVelocity[0]/ pow(r,3);  // setting vxdot = GMx/r^3
    posVel[3] = -4.*M_PI*M_PI*positionAndVelocity[1]/ pow(r,3);  // setting vydot = GMy/r^3

    return posVel;
}   // End of diffEq-function

void energyAndAngMom(vec &positionAndVelocity, ofstream* file)
{
    // Computing velocity squared and radius
    double v = sqrt( pow(positionAndVelocity[2],2) + pow(positionAndVelocity[3],2) );
    double r = sqrt( pow(positionAndVelocity[0],2) + pow(positionAndVelocity[1],2) );

    // Finding energy given velocity and radius
    double EK = v*v/2.;
    double EP = 4*M_PI*M_PI/r;
    double Etot = sqrt(EK*EK + EP*EP);
    double angMom = r*v;

    // Writing energy to file
    *file << setiosflags(ios::showpoint | ios:: uppercase);
    *file << setw(20) << setprecision(15) << EK
          << setw(20) << setprecision(15) << EP
          << setw(20) << setprecision(15) << Etot
          << setw(20) << setprecision(15) << angMom;
}
