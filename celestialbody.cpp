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

void CelestialBody::resetEnergy()
{
    KE = 0;
    PE = 0;
    angMom.setToZero();
}

vec3 CelestialBody::acceleration()
{
    return force/mass;
}

double CelestialBody::totalEnergyOfBody()
{
    return KE + PE;
}
