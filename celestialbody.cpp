#include "celestialbody.h"

CelestialBody::CelestialBody(int index, vec3 position, vec3 velocity, double mass)
{
    this->index     = index;
    this->position  = position;
    this->velocity  = velocity;
    this->mass      = mass;
}

void CelestialBody::resetForce()
{
    force.setToZero();
}
