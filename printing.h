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
    ofstream    positionFile;
    ofstream    velocityFile;
    ofstream    energyAngMomFile;
    string      filenamePrefix;

    // Intitialisation and destructor
    Printing(string filenamePrefix);
    ~Printing();
    void closeAllFiles();

    // Printing functions
    void printingPosition(System &system);
    void printingVelocity(System &system);
    void printingEnergyAngMom(System &system);
    void printingAll(System &system);
    void printingPositionVector(vec3 position);
};


#endif // PRINTING_H
