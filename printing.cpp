#include "printing.h"
#include <vec3.h>

Printing::Printing(string filenamePrefix)
{ // Initialisation of a printer object - it needs a prefix
    this->filenamePrefix = filenamePrefix;
} // End initialisation

Printing::~Printing()
{ // Destructor works by closing all files
    closeAllFiles();
} // End destructor

void Printing::printingPosition(System &system)
{ // Printing position only
    if(!positionFile.is_open()){                                        // Open position file if it's not open
        char *filename = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename, "%s_positions.txt", filenamePrefix.c_str() ); // Create filname with prefix and ending
        positionFile.open(filename);
        delete filename;
    } // End if-statement opening a position file

    for(int i=0; i < system.numberOfBodies(); i++){                     // Print the position of each body
        CelestialBody &body = system.bodies[i];
        positionFile << body.position << " ";
    } // End the for-loop printing the position of each body
    positionFile << std::endl;                                          // Insert a new line when finished
} // End printingPosition-function


void Printing::printingVelocity(System &system)
{ // Printing velocity only
    if(!velocityFile.is_open()) {                                       // Open position file if it's not open
        char *filename = new char[1000];                                // File name can have max 1000 character
        sprintf(filename, "%s_velocities.txt", filenamePrefix.c_str() );// Create filname with prefix and ending
        velocityFile.open(filename);
        delete filename;
    } // End if-statement opening a velocity file

    for(int i=0; i < system.numberOfBodies(); i++){                     // Print the position of each body
        CelestialBody &body = system.bodies[i];
        velocityFile << body.velocity << " ";
    } // End the for-loop printing the position of each body
    velocityFile << std::endl;                                          // Insert a new line when finished
} // End printingVelocity-function


void Printing::printingEnergyAngMom(System &system)
{ // Printing energy and angular momentum only
    if(!energyAngMomFile.is_open()) {                                       // Open position file if it's not open
        char *filename = new char[1000];                                    // File name has max 1000 characters
        sprintf(filename, "%s_energyAngMom.txt", filenamePrefix.c_str() );  // Create filname w prefix and ending
        energyAngMomFile.open(filename);
        delete filename;
    } // End if-statement opening an energy and angular momentum file

    for(int i=0; i < system.numberOfBodies(); i++){                         // Print the position of each body
        CelestialBody &body = system.bodies[i];
        energyAngMomFile << body.velocity << " ";
    } // End the for-loop printing the energy and angular momentum of each body
    energyAngMomFile << std::endl;                                          // Insert a new line when finished
} // End printingEnergyAngMom-function


void Printing::printingAll(System &system)
{ // Function printing position, velocity, energy and angular momentum to file using previously created functions
    printingPosition(system);
    printingVelocity(system);
    printingEnergyAngMom(system);
} // End printingAll-function


void Printing::closeAllFiles()
{ // Function closing all open files
    if(energyAngMomFile.is_open())  energyAngMomFile.close();
    if(positionFile.is_open())      positionFile.close();
    if(velocityFile.is_open())      velocityFile.close();
} // End closeAllFiles-function


void Printing::printingPositionVector(vec3 position)
{ // Printing position only
    if(!positionFile.is_open()){                                        // Open position file if it's not open
        char *filename = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename, "%s_positions.txt", filenamePrefix.c_str() ); // Create filname with prefix and ending
        positionFile.open(filename);
        delete filename;
    } // End if-statement opening a position file

    positionFile << position << " ";
    positionFile << std::endl;                                          // Insert a new line when finished
} // End printingPosition-function
