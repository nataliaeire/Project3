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
        sprintf(filename, "%s_positions.txt", filenamePrefix.c_str() ); // Create filename with prefix and ending
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
        sprintf(filename, "%s_velocities.txt", filenamePrefix.c_str() );// Create filename with prefix and ending
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
    if(!energyAngMomFile.is_open()) {                                           // Open position file if it's not open
        char *filename = new char[1000];                                        // File name has max 1000 characters
        sprintf(filename, "%s_energyAngMom.txt", filenamePrefix.c_str() );      // Create filename w prefix and ending
        energyAngMomFile.open(filename);
        delete filename;
    } // End if-statement opening an energy and angular momentum file

    energyAngMomFile << system.kineticEnergy << " " << system.potentialEnergy << " "
                     << system.totalEnergy() << " " << system.angularMomentum.length() << std::endl;
} // End printingEnergyAngMom-function


void Printing::printingEnergyAngMom(System &system, bool virial)
{ // Printing energy and angular momentum only
    if(virial == true){
        if(!energyAngMomFile.is_open()) {                                       // Open position file if it's not open
            char *filename = new char[1000];                                    // File name has max 1000 characters
            sprintf(filename, "%s_energyAngMom.txt", filenamePrefix.c_str() );  // Create filename w prefix and ending
            energyAngMomFile.open(filename);
            delete filename;
        } // End if-statement opening an energy and angular momentum file

        energyAngMomFile << system.kineticEnergy << " " << system.potentialEnergy << " "
                         << system.totalEnergy() << " " << system.angularMomentum.length() << " "
                         << system.virialKineticEnergy  << " " << system.virialPotentialEnergy << " "
                         << system.boundKineticEnergy+system.boundPotentialEnergy  << " "
                         << system.numberOfBoundBodies() << std::endl;
    }else{
        printingEnergyAngMom(system);
    } // Ending if-statement
} // End printingEnergyAngMom-function


void Printing::printingAll(System &system)
{ // Function printing position, velocity, energy and angular momentum to file using previously created functions
    printingPosition(system);
    printingVelocity(system);
    printingEnergyAngMom(system);
} // End printingAll-function


void Printing::printingPositionXYZ(System &system)
{ // Printing position only
    if(!positionFile.is_open()){                                        // Open position file if it's not open
        char *filename = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename, "%s_positions.xyz", filenamePrefix.c_str() ); // Create filename with prefix and ending
        positionFile.open(filename);
        delete filename;
    } // End if-statement opening a position file

    positionFile << system.numberOfBodies() << std::endl;
    positionFile << "New timestep." << std::endl;
    for(int i=0; i < system.numberOfBodies(); i++){                     // Print the position of each body
        CelestialBody &body = system.bodies[i];
        positionFile << "Ar " << body.position[0] << " " << body.position[1] << " " << body.position[2] << std::endl;
    } // End the for-loop printing the position of each body
} // End printingPosition-function


void Printing::printingPositionXYZ(System &system, int counter)
{ // Printing position only
    if(!positionFile.is_open()){                                        // Open position file if it's not open
        char *filename = new char[1000];                                // File name can have max 1000 characters
        sprintf(filename, "%s_positions.xyz", filenamePrefix.c_str() ); // Create filename with prefix and ending
        positionFile.open(filename);
        delete filename;
    } // End if-statement opening a position file

    positionFile << system.numberOfBodies() << std::endl;
    positionFile << "Time step " << counter << "." << std::endl;
    for(int i=0; i < system.numberOfBodies(); i++){                     // Print the position of each body
        CelestialBody &body = system.bodies[i];
        positionFile << "Ar " << body.position[0] << " " << body.position[1] << " " << body.position[2] << std::endl;
    } // End the for-loop printing the position of each body
} // End printingPosition-function


void Printing::printingAll(System &system, int counter, int n)
{ // Function printing only each n'th position, velocity, energy and angular momentum to file using previously created functions
    if (counter % n == 0)           printingAll(system);
} // End printingAll-function


void Printing::closeAllFiles()
{ // Function closing all open files
    if(energyAngMomFile.is_open())  energyAngMomFile.close();
    if(positionFile.is_open())      positionFile.close();
    if(velocityFile.is_open())      velocityFile.close();
} // End closeAllFiles-function


void Printing::printing3Vector(vec3 vector, std::string filenameEnding)
{ // Printing position only
    this->filenameEnding = filenameEnding;
    if(!vectorFile.is_open()){                                                              // Open position file if it's not open
        char *filename = new char[1000];                                                    // File name can have max 1000 characters
        sprintf(filename, "%s_%s.txt", filenamePrefix.c_str(), filenameEnding.c_str() );    // Create filename with prefix and ending
        vectorFile.open(filename);
        delete filename;
    } // End if-statement opening a position file
    vectorFile << vector << " " << std::endl;
} // End printingPosition-function
