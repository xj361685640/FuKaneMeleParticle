#include "Shell.h"

Shell::Shell (double pSize, double orderpSize) {
    radius = pSize; // radius in distance units
    innerRadius = orderpSize; // used to specify inner radius of the shell
    shape = "Shell";
}

bool Shell::WithinParticle (rowvec dummyVec) {
    if (dot(dummyVec,dummyVec) < radius*radius && dot(dummyVec,dummyVec) > innerRadius*innerRadius) {
        return true;
    } else {
        return false;
    }
}

void Shell::AddDisorder (double disorderStrength, double disorderCoverage, double delta){ 
	cout << "Disorder not implemented for the shell" << endl;
} 

void Shell::PrintInfo (double disorderStrength) {
    cout << "Total # of sublattice sites: " << 2*cells << " within shell of inner radius " << innerRadius  << " and outer radius " << radius << endl;
}
