#include "printing.h"
#include <vec3.h>

Printing::Printing(string filenamePrefix)
{
    this->filenamePrefix = filenamePrefix;
}

Printing::~Printing()
{
    closeAllFiles();
}

void Printing::printingPosition(System &system)
{
    if(!positionFile.is_open()) {
        char *filename = new char[1000];
        sprintf(filename, "%s_positions.txt", filenamePrefix.c_str() );
        positionFile.open(filename);
        delete filename;
    }

    for(int i=0; i < system.numberOfBodies(); i++){
        CelestialBody &body = system.bodies[i];
        positionFile << body.position << " ";
    }
    positionFile << std::endl;
}


void Printing::printingVelocity(System &system)
{
    if(!velocityFile.is_open()) {
        char *filename = new char[1000];
        sprintf(filename, "%s_velocities.txt", filenamePrefix.c_str() );
        velocityFile.open(filename);
        delete filename;
    }

    for(int i=0; i < system.numberOfBodies(); i++){
        CelestialBody &body = system.bodies[i];
        velocityFile << body.velocity << " ";
    }
    velocityFile << std::endl;
}


void Printing::printingEnergyAngMom(System &system)
{
    if(!energyAngMomFile.is_open()) {
        char *filename = new char[1000];
        sprintf(filename, "%s_energyAngMom.txt", filenamePrefix.c_str() );
        energyAngMomFile.open(filename);
        delete filename;
    }

    for(int i=0; i < system.numberOfBodies(); i++){
        CelestialBody &body = system.bodies[i];
        energyAngMomFile << body.velocity << " ";
    }
    energyAngMomFile << std::endl;
}


void Printing::printingAll(System &system)
{
    printingPosition(system);
    printingVelocity(system);
    printingEnergyAngMom(system);
}


void Printing::closeAllFiles()
{
    if(energyAngMomFile.is_open())  energyAngMomFile.close();
    if(positionFile.is_open())      positionFile.close();
    if(velocityFile.is_open())      velocityFile.close();
}
