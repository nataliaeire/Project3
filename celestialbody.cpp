#include "celestialbody.h"

CelestialBody::CelestialBody(int index, vec3 position, vec3 velocity, double mass)
{  //  Initialising the body
    this->index     = index;
    this->position  = position;
    this->velocity  = velocity;
    this->mass      = mass;
    gravitationallyBound = true;
}

void CelestialBody::resetForce()
{ // Setting the force on the body to zero
    force.setToZero();
}

void CelestialBody::resetEnergy()
{  // Resetting the kinetic, potential ang angular momentum of a body
    KE = 0;
    PE = 0;
    angMom.setToZero();
}

vec3 CelestialBody::acceleration()
{   // Setting the acceleration
    return force/mass;
}

double CelestialBody::totalEnergyOfBody()
{  // Finding the total energy of the body for use in the cluser problem
    return KE + PE;
}
