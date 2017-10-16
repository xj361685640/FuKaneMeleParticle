#ifndef GUARD_SPHERICAL_H
#define GUARD_SPHERICAL_H

#include <iostream>
#include "Particle.h"
class Spherical: public Particle {
    public:
        double radius;
        double orderRadius;
        Spherical (double pSize, double orderpSize);

        bool WithinParticle (rowvec dummyVec);

        void AddDisorder (double disorderStrength, double disorderCoverage, double delta);

        void PrintInfo (double disorderStrength);
};

#endif
