#ifndef GUARD_SHELL_H
#define GUARD_SHELL_H

#include <iostream>
#include "Particle.h"
class Shell: public Particle {
    public:
		double innerRadius;
        double radius;
        Shell (double pSize, double orderpSize);

        bool WithinParticle (rowvec dummyVec);

        void AddDisorder (double disorderStrength, double disorderCoverage, double delta);

        void PrintInfo (double disorderStrength);
};

#endif
