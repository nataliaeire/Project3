#ifndef PRINTING_H
#define PRINTING_H
#include <fstream>
#include <vector>
#include <string>
#include <vec3.h>
#include <celestialbody.h>
#include <system.h>
using std::ofstream; using std::string;

class Printing
{
public:
    ofstream    vectorFile;
    ofstream    positionFile;
    ofstream    positionxyzFile;
    ofstream    velocityFile;
    ofstream    energyAngMomFile;
    string      filenamePrefix;
    string      filenameEnding;

    // Intitialisation and destructor
    Printing(string filenamePrefix);
    ~Printing();
    void closeAllFiles();

    // Printing functions
    void printingPosition(System &system);
    void printingPosition(System &system, bool virial);
    void printingVelocity(System &system);
    void printingEnergyAngMom(System &system);
    void printingEnergyAngMom(System &system, bool virial);
    void printingAll(System &system);
    void printingAll(System &system, int counter, int n);
    void printingAll(System &system, int counter, bool virial);
    void printing3Vector(vec3 vector, string filenameEnding);
    void printingPositionXYZ(System &system);
    void printingPositionXYZ(System &system, int counter);
};   // End of Printing class declaration


#endif // PRINTING_H
