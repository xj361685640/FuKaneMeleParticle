#ifndef GUARD_RHOMBOHEDRAL_H
#define GUARD_RHOMBOHEDRAL_H

#include <iostream>
#include "Particle.h"

class Rhombohedral: public Particle {
    public:
        double edgeSize;
        double orderEdgeSize;
        Rhombohedral (double pSize, double orderpSize);

        bool WithinParticle (rowvec dummyVec);

        void PrintInfo (double disorderStrength);
        
        void AddDisorder(double disorderStrength, double disorderCoverage, double delta);
};

#endif
