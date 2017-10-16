#include "Spherical.h"

Spherical::Spherical (double pSize, double orderpSize) {
    radius = pSize; // radius in distance units
    orderRadius = orderpSize; // radius above which disorder starts
    shape = "Spherical";
}

bool Spherical::WithinParticle (rowvec dummyVec) {
    if (dot(dummyVec,dummyVec) < radius*radius) {
        return true;
    } else {
        return false;
    }
}

void Spherical::AddDisorder (double disorderStrength, double disorderCoverage, double delta){ 
    double dummy;
    rowvec dummyVec(3);
    disorderedSites=0;
    srand(time(0));
    for (int i=0; i<cells; i++) {
        if (1.0*rand()/RAND_MAX < disorderCoverage) {
            dummyVec=sublatOne(i, span(1,3));
            if (dot(dummyVec,dummyVec)>orderRadius*orderRadius-delta) {
                disorderedSites++;
                dummy=(-1.0+2.0*rand()/RAND_MAX);
                Hamiltonian(4*i,4*i)=disorderStrength*dummy;
                Hamiltonian(4*i+1,4*i+1)=disorderStrength*dummy;
            }
            dummyVec=sublatTwo(i, span(1,3));
            if (dot(dummyVec,dummyVec)>orderRadius*orderRadius-delta) {
                disorderedSites++;
                dummy=(-1.0+2.0*rand()/RAND_MAX);
                Hamiltonian(4*i+2,4*i+2)=disorderStrength*dummy;
                Hamiltonian(4*i+3,4*i+3)=disorderStrength*dummy;
            }
        }
    }
} 

void Spherical::PrintInfo (double disorderStrength) {
    cout << "Total # of sublattice sites: " << 2*cells << " within sphere of radius " << radius  << endl;
    cout << "# of disordered sublattice sites: " << disorderedSites << " Strength: " << disorderStrength << " above radius " << orderRadius << endl;
}
