#include "Rhombohedral.h"

Rhombohedral::Rhombohedral (double pSize, double orderpSize) {
    edgeSize = pSize;
    orderEdgeSize = orderpSize;
    shape = "Rhombohedral";
}

bool Rhombohedral::WithinParticle (rowvec dummyVec) {
    return true;
} 

void Rhombohedral::PrintInfo (double disorderStrength) {
    cout << "Total # of sublattice sites: " << 2*cells << " within rhomboid of size " << edgeSize << endl;
    cout << "# of disordered sublattice sites: " << disorderedSites << " Strength: " << disorderStrength << " coveing the whole rhombohedral particle" << endl;
}

void Rhombohedral::AddDisorder (double disorderStrength, double disorderCoverage, double delta){ 
    double dummy;
    disorderedSites=0;
    srand(time(0));
    for (int i=0; i<cells; i++) {
        if (1.0*rand()/RAND_MAX < disorderCoverage) {
            disorderedSites++;
            dummy=(-1.0+2.0*rand()/RAND_MAX);
            Hamiltonian(4*i,4*i)=disorderStrength*dummy;
            Hamiltonian(4*i+1,4*i+1)=disorderStrength*dummy;
            disorderedSites++;
            dummy=(-1.0+2.0*rand()/RAND_MAX);
            Hamiltonian(4*i+2,4*i+2)=disorderStrength*dummy;
            Hamiltonian(4*i+3,4*i+3)=disorderStrength*dummy;
        }
    }
}

