#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "KeplerClasses.h"

using namespace std;
using namespace arma;


int main(){

    //CelestialBody earth(1, 0, 0, 2*M_PI, 0, 0, 3e-6);
    System solarsyst;
    vec pos(2);
    pos[0] = 1;
    pos[1] = 0;
    vec vel(2);
    vel[0] = 0;
    vel[1] = 2*M_PI;
    double mass = 3e-6;

    solarsyst.addBody(pos,vel,mass);
    solarsyst.definingBodies();
//    vec ourBody = solarsyst.defineSystem;
//    cout << ourBody << endl;

    return 0;
}
