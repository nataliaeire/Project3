#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace arma;

vec diffEq(vec positionAndVelocity, vec masses, int numberOfBodies, ofstream *enmom);
void rk4(double h, int numberOfBodies, vec masses, vec &positionAndVelocity, ofstream *posvel, ofstream *enmom);
void calculateForcesAndEnergy(int numberOfBodies, vec positionAndVelocity,
                              vec &energyAngMom, vec &Fx, vec &Fy, vec masses);
//void energyAndAngMom(vec &positionAndVelocity, ofstream* file);

int main()
{
    int n = 1000;                  // Number of iterations
    double T = 1;                  // Total time of simulation (yrs)
    double h = T/n;                // Step size
    int numberOfBodies = 2;

    // Create vector containting all positions and velocities
    vec positionAndVelocity(4*numberOfBodies);
    positionAndVelocity(0) = 1;     // Initial x-position
    positionAndVelocity(1) = 0;     // Initial y-position
    positionAndVelocity(2) = 0;     // Initial vx-velocity
    positionAndVelocity(3) = 2*M_PI;  // Initial vy-velocity (2*1 AU*pi / 1 yr)
    positionAndVelocity(4) = 0;     // Initial x-position
    positionAndVelocity(5) = 0;     // Initial y-position
    positionAndVelocity(6) = 0;     // Initial vx-velocity
    positionAndVelocity(7) = 0;     // Initial vy-velocity (2*1 AU*pi / 1 yr)

    vec masses(numberOfBodies);
    masses(0) = 5.972*pow(10,24)/(1.9891*pow(10,30));
    masses(1) = 1.;

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
        //energyAndAngMom(positionAndVelocity, &energyOut);
        rk4(h, numberOfBodies, masses, positionAndVelocity, &posVelOut, &energyOut);
    }

    // Closing files
    posVelOut.close();
    energyOut.close();

    return 0;
}   // End of main function


// Function performing the Runge-Kutta 4 method
void rk4(double h, int numberOfBodies, vec masses, vec &positionAndVelocity, ofstream* posvel, ofstream* enmom)
{
    // Declaring slopes to perform RK4
    vec k1, k2, k3, k4;

    k1 = h*diffEq(positionAndVelocity,       masses, numberOfBodies, enmom);
    k2 = h*diffEq(positionAndVelocity+k1/2., masses, numberOfBodies, enmom);
    k3 = h*diffEq(positionAndVelocity+k2/2., masses, numberOfBodies, enmom);
    k4 = h*diffEq(positionAndVelocity+k3,    masses, numberOfBodies, enmom);

    // Updating position and velocity after RK4
    positionAndVelocity = positionAndVelocity + (k1 + 2.*k2 + 2.*k3 + k4)/6.;

    // Writing position and velocity to file
    *posvel << setiosflags(ios::showpoint | ios:: uppercase);
    *posvel << setw(20) << setprecision(15) << positionAndVelocity.t();
}


// Function setting the differential equations
vec diffEq(vec positionAndVelocity, vec masses, int numberOfBodies, ofstream* enmom)
{
    vec dotPosVel(4*numberOfBodies);                // Creating temporary position and velocity vector
    vec Fx = zeros(numberOfBodies);                 // Gravitational force in x-direction to compute vxdot
    vec Fy = zeros(numberOfBodies);                 // Gravitational force in y-direction to compute vydot
    vec energyAngMom = zeros(4);
    calculateForcesAndEnergy(numberOfBodies, positionAndVelocity, energyAngMom, Fx, Fy, masses);

    // Writing energy to file
    *enmom << setiosflags(ios::showpoint | ios:: uppercase);
    *enmom << setw(20) << setprecision(15) << energyAngMom.t();

    // Finding the derivative of each velocity from the gravitational force
    for(int i = 0; i < numberOfBodies; i++){
        dotPosVel[4*i+0] = positionAndVelocity[4*i+2];  // setting xdot = vx
        dotPosVel[4*i+1] = positionAndVelocity[4*i+3];  // setting ydot = vy
        dotPosVel[4*i+2] = Fx[i]/masses[i];
        dotPosVel[4*i+3] = Fy[i]/masses[i];
    }   // Ending for-loop computing dotv

    return dotPosVel;
}   // End of diffEq-function


// Function calculating forces and energy (and angular momentum!) for the system
void calculateForcesAndEnergy(int numberOfBodies, vec positionAndVelocity,
                              vec &energyAngMom, vec &Fx, vec &Fy, vec masses)
{
    double G = 4*M_PI*M_PI;                         // Defining the gravitational constant
    energyAngMom.zeros();
    Fx.zeros();
    Fy.zeros();
    for(int i = 0; i < numberOfBodies; i++)
    {
        for(int j=i+1; j < numberOfBodies; j++){
            double dx = positionAndVelocity[4*i+0]-positionAndVelocity[4*j+0];
            double dy = positionAndVelocity[4*i+1]-positionAndVelocity[4*j+1];
            double dr2 = dx*dx + dy*dy;
            double dr = sqrt(dr2);

            // Updating gravitational force experienced by celestial object
            double factor = G*masses(j)*masses(i) / pow(dr,3);
            energyAngMom(1) += factor*dr;
            Fx[i] -= factor*dx;
            Fy[i] -= factor*dy;
            Fx[j] += factor*dx;
            Fy[j] += factor*dy;
        }   // Ending for-loop computing force

        double v2 = pow(positionAndVelocity[4*i+2],2) + pow(positionAndVelocity[4*i+3],2);
        double pos = sqrt(pow(positionAndVelocity[4*i+0],2) + pow(positionAndVelocity[4*i+1],2));
        energyAngMom(0) += 0.5*masses(i)*v2;
        energyAngMom(3) += masses(i)*pos*sqrt(v2);
        energyAngMom(2) = sqrt(pow(energyAngMom(0),2) + pow(energyAngMom(1),2));

    }   // Ending for-loop going over all celestial bodies
}   // Ending calculateForcesAndEnergy-function

