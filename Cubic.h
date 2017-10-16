#ifndef GUARD_CUBIC_H
#define GUARD_CUBIC_H

#include <iostream>
#include "Particle.h"
class Cubic: public Particle {
    public: 
        double edgeSize;
        double orderEdgeSize;
        Cubic (double pSize, double orderpSize);

        bool WithinParticle (rowvec dummyVec); 

        void AddDisorder(double disorderStrength, double disorderCoverage, double delta);
    
        void PrintInfo (double disorderStrength);
};
#endif
