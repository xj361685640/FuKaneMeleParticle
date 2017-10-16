#include "Cubic.h"

Cubic::Cubic (double pSize, double orderpSize) {
    edgeSize = pSize;
    orderEdgeSize = orderpSize;
    shape = "Cubic";
}

bool Cubic::WithinParticle (rowvec dummyVec) {
    if (dummyVec(0) > -edgeSize && dummyVec(0) < edgeSize &&
            dummyVec(1) > -edgeSize && dummyVec(1) < edgeSize &&
            dummyVec(2) > -edgeSize && dummyVec(2) < edgeSize ) {
        return true;
    } else {
        return false;
    }
}

void Cubic::AddDisorder(double disorderStrength, double disorderCoverage, double delta){ 
    double dummy;
    rowvec dummyVec(3);
    disorderedSites=0;
    srand(time(0));
    for (int i=0; i<cells; i++) {
        if (1.0*rand()/RAND_MAX < disorderCoverage) {
            dummyVec=sublatOne(i, span(1,3)); 
            if (dummyVec(0) > -orderEdgeSize && dummyVec(0) < orderEdgeSize &&
                    dummyVec(1) > -orderEdgeSize && dummyVec(1) < orderEdgeSize &&
                    dummyVec(2) > -orderEdgeSize && dummyVec(2) < orderEdgeSize ) {
                disorderedSites++;
                dummy=(-1.0+2.0*rand()/RAND_MAX);
                Hamiltonian(4*i,4*i)=disorderStrength*dummy;
                Hamiltonian(4*i+1,4*i+1)=disorderStrength*dummy;
            }
            dummyVec=sublatTwo(i, span(1,3)); 
            if (dummyVec(0) > -orderEdgeSize && dummyVec(0) < orderEdgeSize &&
                    dummyVec(1) > -orderEdgeSize && dummyVec(1) < orderEdgeSize &&
                    dummyVec(2) > -orderEdgeSize && dummyVec(2) < orderEdgeSize ) {
                disorderedSites++;
                dummy=(-1.0+2.0*rand()/RAND_MAX);
                Hamiltonian(4*i+2,4*i+2)=disorderStrength*dummy;
                Hamiltonian(4*i+3,4*i+3)=disorderStrength*dummy;
            }
        }
    }
} 

void Cubic::PrintInfo (double disorderStrength) {
    cout << "Total # of sublattice sites: " << 2*cells << " within cube of size " << edgeSize << endl;
    cout << "# of disordered sublattice sites: " << disorderedSites << " Strength: " << disorderStrength << " in cubic shell of size > " << orderEdgeSize << endl;
}



