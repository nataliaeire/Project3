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
    void printingVelocity(System &system);
    void printingEnergyAngMom(System &system);
    void printingAll(System &system);
    void printing3Vector(vec3 vector, string filenameEnding);
};


#endif // PRINTING_H
