#include "celestialbody.h"
using namespace arma;

CelestialBody::CelestialBody(vec position, vec velocity, double mass)
{
    this->position    = position;
    this->velocity    = velocity;
    this->mass        = mass;

    this->positionDot = zeros(3);
    this->velocityDot = zeros(3);
    this->force       = zeros(3);

    resetForce();
}

void CelestialBody::resetForce()
{
    force.zeros();
}
