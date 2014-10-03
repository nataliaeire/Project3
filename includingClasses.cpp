#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace arma;


// ================================= CLASS: CelesitalBody ============================================== //
class CelestialBody {
public:
    double x;
    double y;
    double vx;
    double vy;
    double mass;

    CelestialBody();
    CelestialBody (double x, double y, double vx, double vy, double mass);
};

CelestialBody::CelestialBody () {
}

CelestialBody::CelestialBody (double x, double y, double vx, double vy, double mass){
    this->x = x;
    this->y = y;
    this->vx = vx;
    this->vy = vy;
    this->mass = mass;
}


// ================================= CLASS: SYSTEM ============================================== //
class System {
public:
    int nBodies;
    CelestialBody* bodies;

    System(int n);
    vec defineBody;

    void twoBodies();
    void threeBodies();
    void definingBodies();

    void RK4(double h, int nBodies, ofstream* posvel, ofstream* enmom);
    vec  diffEq(vec, int nBodies, ofstream*);
};

System::System(int n) {
    this->nBodies = n;
    this->bodies = new CelestialBody[nBodies];

    for (int i = 0; i < nBodies; i ++) {
        this->bodies[i] = CelestialBody ();
    }
}

void System::twoBodies() {
    // Sola
    this->bodies[0].mass = 1.0;
    this->bodies[0].x = 0.0;
    this->bodies[0].y = 0.0;
    this->bodies[0].vx = 0.0;
    this->bodies[0].vy = 0.0;

    // Jorda
    this->bodies[1].mass = 3e-6;
    this->bodies[1].x = 1.0;
    this->bodies[1].y = 0.0;
    this->bodies[1].vx = 2*M_PI;
    this->bodies[1].vy = 0.0;
}

void System::threeBodies() {
    // Sola
    this->bodies[0].mass = 1.0;
    this->bodies[0].x = 0.0;
    this->bodies[0].y = 0.0;
    this->bodies[0].vx = 0.0;
    this->bodies[0].vy = -M_PI*6e-6;

    // Jorda
    this->bodies[1].mass = 3e-6;
    this->bodies[1].x = 1.0;
    this->bodies[1].y = 0.0;
    this->bodies[1].vx = 0.0;
    this->bodies[1].vy = 2*M_PI;

    // Jupiter
    this->bodies[2].mass = 1e-3;
    this->bodies[2].x = 0.0;
    this->bodies[2].y = 5.20;
    this->bodies[2].vx = 0.0;
    this->bodies[2].vy = 0.0;
}

void System::definingBodies() {
    defineBody(5*2); //this should be 5*nBodies

    for(int i = 0; i < 2; i+=5) {
        defineBody[5*i]   = bodies[i].x;
        defineBody[5*i+1] = bodies[i].y;
        defineBody[5*i+2] = bodies[i].vx;
        defineBody[5*i+3] = bodies[i].vy;
        defineBody[5*i+4] = bodies[i].mass;
    }
}


/*
void System::nBodies(int n) {
    for (int i=1..n)
        this->bodies[i].x = rand();
}*/

/*
// Function setting the differential equations
vec System::diffEq(vec k, ofstream* enmom)
{
    vec dotPosVel(5*nBodies);                // Creating temporary position and velocity vector
    vec Fx = zeros(nBodies);                 // Gravitational force in x-direction to compute vxdot
    vec Fy = zeros(nBodies);                 // Gravitational force in y-direction to compute vydot
    vec energyAngMom = zeros(5);
    calculateForcesAndEnergy(nBodies, positionAndVelocity, energyAngMom, Fx, Fy, masses);

    // Writing energy to file
    *enmom << setiosflags(ios::showpoint | ios:: uppercase);
    *enmom << setw(20) << setprecision(15) << energyAngMom.t();

    // Finding the derivative of each velocity from the gravitational force
    for(int i = 0; i < nBodies; i++){

        //posVel;

        dotPosVel[5*i+0] = positionAndVelocity[4*i+2];  // setting xdot = vx
        dotPosVel[5*i+1] = positionAndVelocity[4*i+3];  // setting ydot = vy
        dotPosVel[5*i+2] = Fx[i]/masses[i];
        dotPosVel[5*i+3] = Fy[i]/masses[i];
    }   // Ending for-loop computing dotv

    return dotPosVel;
}   // End of diffEq-function
*/

vec System::diffEq(vec, int, ofstream*) {
    vec hei = vec(3);
    return hei;
}

// Function performing the Runge-Kutta 4 method
void System::RK4(double h, int nBodies, ofstream* posvel, ofstream* enmom)
{
    // Declaring slopes to perform RK4
    vec k0, k1, k2, k3, k4;

    k0 = zeros(5*nBodies);
    k1 = h*diffEq(k0,    nBodies, enmom);
    k2 = h*diffEq(k1/2., nBodies, enmom);
    k3 = h*diffEq(k2/2., nBodies, enmom);
    k4 = h*diffEq(k3,    nBodies, enmom);

    // Updating position and velocity after RK4
    defineBody = defineBody + (k1 + 2.*k2 + 2.*k3 + k4)/6.;

    // Writing position and velocity to file
    *posvel << setiosflags(ios::showpoint | ios:: uppercase);
    *posvel << setw(20) << setprecision(15) << defineBody.t();
}




//// Function calculating forces and energy (and angular momentum!) for the system
//void calculateForcesAndEnergy(int nBodies, vec positionAndVelocity,
//                              vec &energyAngMom, vec &Fx, vec &Fy, vec masses)
//{
//    double G = 4*M_PI*M_PI;                         // Defining the gravitational constant
//    energyAngMom.zeros();
//    Fx.zeros();
//    Fy.zeros();
//    for(int i = 0; i < nBodies; i++)
//    {
//        for(int j=i+1; j < nBodies; j++){
//            double dx = positionAndVelocity[4*i+0]-positionAndVelocity[4*j+0];
//            double dy = positionAndVelocity[4*i+1]-positionAndVelocity[4*j+1];
//            double dr2 = dx*dx + dy*dy;
//            double dr = sqrt(dr2);

//            // Updating gravitational force experienced by celestial object
//            double factor = G*masses(j)*masses(i) / pow(dr,3);
//            energyAngMom(1) += factor*dr;
//            Fx[i] -= factor*dx;
//            Fy[i] -= factor*dy;
//            Fx[j] += factor*dx;
//            Fy[j] += factor*dy;
//        }   // Ending for-loop computing force

//        double v2 = pow(positionAndVelocity[4*i+2],2) + pow(positionAndVelocity[4*i+3],2);
//        double pos = sqrt(pow(positionAndVelocity[4*i+0],2) + pow(positionAndVelocity[4*i+1],2));
//        energyAngMom(0) += 0.5*masses(i)*v2;
//        energyAngMom(3) += masses(i)*pos*sqrt(v2);
//        energyAngMom(2) = sqrt(pow(energyAngMom(0),2) + pow(energyAngMom(1),2));

//    }   // Ending for-loop going over all celestial bodies
//}   // Ending calculateForcesAndEnergy-function









// ================================= MAIN FUNCTION ============================================== //
int mainInClasses(){

    double T = 50;                  // Total time of simulation (yrs)
    double h = 1/12;                // Step size
    int n = T/h;                    // Number of iterations
    int nBodies = 2;


    return 0;
}
