#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace arma;

vec diffEq(vec positionAndVelocity);
void rk4(double h, vec &positionAndVelocity, ofstream *file);
void energyAndAngMom(vec &positionAndVelocity, ofstream* file);

int main()
{
    int n = 1000;                   // Number of iterations
    double T = 1;                  // Total time of simulation (yrs)
    double h = T/n;                 // Step size

    // Create vector containting all positions and velocities
    vec positionAndVelocity(4);
    positionAndVelocity(0) = 1;     // Initial x-position
    positionAndVelocity(1) = 0;     // Initial y-position
    positionAndVelocity(2) = 0;     // Initial vx-velocity
    positionAndVelocity(3) = 2*M_PI;  // Initial vy-velocity (2*1 AU*pi / 1 yr)

    // Creating file to write to
    ofstream posVelOut;
    ofstream energyOut;
    posVelOut.open("posvel.dat");
    energyOut.open("energyangmom.dat");
    // Writing initial positions and velocities to file
    posVelOut << setiosflags(ios::showpoint | ios:: uppercase);
    posVelOut << setw(20) << setprecision(15) << positionAndVelocity.t();


    // Running RK4
    for(int i=0; i<n; i++){
        energyAndAngMom(positionAndVelocity, &energyOut);
        rk4(h, positionAndVelocity, &posVelOut);
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

// Function setting the differential equations
vec diffEq(vec positionAndVelocity)
{
    vec posVel(4);      // Creating temporary position and velocity vector
    double r = sqrt( pow(positionAndVelocity[0],2)+pow(positionAndVelocity[1],2) );   // Calculating radius

    posVel[0] = positionAndVelocity[2];                         // setting xdot = vx
    posVel[1] = positionAndVelocity[3];                         // setting ydot = vy
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
