#include "celestialbody.h"

CelestialBody::CelestialBody(vec3 position, vec3 velocity, double mass)
{
    this->position  = position;
    this->velocity  = velocity;
    this->mass      = mass;
}

void CelestialBody::resetForce()
{
    force.setToZero();
}
