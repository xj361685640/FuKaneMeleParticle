#include "Particle.h"

Particle::Particle() {}
Particle::~Particle() {}   

void Particle::SetDataStructures (int pCells) {
    cells = pCells;
    sublatOne = zeros(cells,4);
    sublatTwo = zeros(cells,4);
    Hamiltonian.set_size(4*cells, 4*cells);
    Hamiltonian.set_real(zeros(4*cells, 4*cells));
    Hamiltonian.set_imag(zeros(4*cells, 4*cells));
    eigvecs.set_size(4*cells, 4*cells);
    eigvecs.set_real(zeros(4*cells, 4*cells));
    eigvecs.set_imag(zeros(4*cells, 4*cells));
    eigvals = zeros(4*cells);
}

bool Particle::WithinParticle (rowvec dummyVec) {
    return 0;
}

void Particle::AddDisorder(double disorderStrength, double disorderCoverage, double delta) {}

void Particle::PrintInfo(double disorderStrength) {}
