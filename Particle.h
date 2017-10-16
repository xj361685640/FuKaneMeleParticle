#ifndef GUARD_PARTICLE_H
#define GUARD_PARTICLE_H
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

class Particle {
    public:  
        Particle();
        ~Particle();

        int cells;
        int disorderedSites;
        string shape;
        mat sublatOne, sublatTwo;
        cx_mat Hamiltonian;
        cx_mat eigvecs;
        vec eigvals;

        void SetDataStructures (int);
        virtual bool WithinParticle (rowvec);
        virtual void AddDisorder(double, double, double);
        virtual void PrintInfo(double);
};

#endif
