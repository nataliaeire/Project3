#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H
#include <armadillo>

class CelestialBody{
public:
    arma::vec position;
    arma::vec velocity;
    arma::vec force;
    double mass;

    CelestialBody(arma::vec position, arma::vec velocity, double mass);
    void resetForce();

};  // End of CelestialBody class declaration


#endif // CELESTIALBODY_H
