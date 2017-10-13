#include <iostream>
#include <armadillo>                                                                                 
#include<complex>                                                                                    
#include<cstdlib>
#include <ctime>
#include <string>
using namespace std;                                                                                 
using namespace arma;                                                                                

class Particle {
    public:
        Particle() {}
        ~Particle() {}   

        int cells;
        int disorderedSites;
        string shape;
        mat sublatOne, sublatTwo;
        cx_mat Hamiltonian;
        cx_mat eigvecs;
        vec eigvals;

        void SetDataStructures (int pCells) {
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
};

class Spherical: public Particle {
    private:
        double dummy;
    public:
        double radius;
        double orderRadius;
        Spherical(double pSize, double orderpSize) {
            radius = pSize; // radius in distance units
            orderRadius = orderpSize; // radius above which disorder starts
            shape = "Spherical";
        }

        bool WithinParticle (rowvec dummyVec) {
            if (dot(dummyVec,dummyVec) < radius*radius) {
                return true;
            } else {
                return false;
            }
        }

        bool OutsideOrderCore (rowvec dummyVec) {
            if (dot(dummyVec, dummyVec) > orderRadius*orderRadius) {
                return true;
            } else {
                return false;
            }
        }

        void AddDisorder(double disorderStrength, double disorderCoverage, double delta){ 
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

        void PrintInfo (double disorderStrength) {
            cout << "Total # of sublattice sites: " << 2*cells << " within radius " << radius << " of " << shape << "particle" << endl;
            cout << "# of disordered sublattice sites: " << disorderedSites << " Strength: " << disorderStrength << " above radius " << orderRadius << endl;
        }
};

class Cubic: public Particle {
    public: 
        double edgeSize;
        double orderEdgeSize;
        Cubic (double pSize, double orderpSize) {
            edgeSize = pSize;
            orderEdgeSize = orderpSize;
            shape = "Cubic";
        }

        bool WithinParticle (rowvec dummyVec) {
            if (dummyVec(0) > -edgeSize && dummyVec(0) < edgeSize &&
                    dummyVec(1) > -edgeSize && dummyVec(1) < edgeSize &&
                    dummyVec(2) > -edgeSize && dummyVec(2) < edgeSize ) {
                return true;
            } else {
                return false;
            }
        }

        bool OutsideOrderCore (rowvec dummyVec) {
            if (dummyVec(0) < -orderEdgeSize && dummyVec(0) > orderEdgeSize &&
                    dummyVec(1) < -orderEdgeSize && dummyVec(1) > orderEdgeSize &&
                    dummyVec(2) < -orderEdgeSize && dummyVec(2) > orderEdgeSize ) {
                return true;
            } else {
                return false;
            }
        }

        void PrintInfo (double disorderStrength) {
            cout << "Total # of sublattice sites: " << 2*cells << " within cube of size " << edgeSize << endl;
            cout << "# of disordered sublattice sites: " << disorderedSites << " Strength: " << disorderStrength << " in cubic shell of size > " << orderEdgeSize << endl;
        }

};

class Rhombohedral: public Particle {
    public:
        double edgeSize;
        double orderEdgeSize;
        Rhombohedral (double pSize, double orderpSize) {
            edgeSize = pSize;
            orderEdgeSize = orderpSize;
            shape = "Rhombohedral";
        }

        bool WithinParticle (rowvec dummyVec) {
            return true;
        } 
        void PrintInfo (double disorderStrength) {
            cout << "Total # of sublattice sites: " << 2*cells << " within cube of size " << edgeSize << endl;
            cout << "# of disordered sublattice sites: " << disorderedSites << " Strength: " << disorderStrength << " in cubic shell of size > " << orderEdgeSize << endl;
        }

        bool OutsideOrderCore (rowvec dummyVec) { 
            return true; // TODO
        }


};


int main(){
    // Input 
    int sizeCells;  //sizeCells of lattice in +/-ve x/y/z direction
    int pSizeCells;  // size of the particle (radius, edge) in unit cells - carved out of bulk lattice
    double disorderStrength; // Magnitude of disorder potential
    double disorderCoverage; // fraction of disordered states above orderRadius

    // Parameters and data structures (Don't need to change those)
    const double tHopping = 1;
    const double deltatHopping = 0.4; 
    const complex<double> spinOrbitCoupling(0,1); 

    // Geometry of the model;
    const rowvec latVec1 ("0 0.5 0.5");
    const rowvec latVec2 ("0.5 0 0.5");
    const rowvec latVec3 ("0.5 0.5 0");
    const rowvec sublatVec ("0.25 0.25 0.25");
    const rowvec comShift ("-0.5 -0.5 -0.5"); // to have centre of mass of sublatOne at the origin
    const double latVecNorm = 1/sqrt(2); // cubic cell size a = 1
    const double sublatVecNorm = sqrt(3)/4;

    const double delta = 0.01; // small number for numerical comparisons

    // TODO do input file with keywords
    ifstream inputfile;
    inputfile.open("input.txt");
    while (inputfile >> sizeCells >> pSizeCells >> disorderStrength >> disorderCoverage) 
        cout << "Input: "<< sizeCells << " " << pSizeCells << " " <<  disorderStrength << " " << disorderCoverage << endl;

    //Spherical myParticle(pSizeCells*latVecNorm, pSizeCells*latVecNorm - sublatVecNorm);
    //Cubic myParticle(pSizeCells*latVecNorm, pSizeCells*latVecNorm - sublatVecNorm);
    Rhombohedral myParticle(pSizeCells*latVecNorm, pSizeCells*latVecNorm - sublatVecNorm);

    int nCells = 8*sizeCells*sizeCells*sizeCells; // # unit cells 
    int pCells = 0; // # unit cells belonging to the particle

    int dummy = 0, index = 0; 
    rowvec dummyVec(3);
    cx_mat dummyMat(2,2);
    rowvec pcomShift(3);

    mat sublatOne(nCells, 4); // Sites at origin of a unit cell
    mat sublatTwo(nCells, 4); // Sites at sublatVec in a unit cell

    cx_mat unity(2,2), sigma_x(2,2), sigma_y(2,2), sigma_z(2,2); // Pauli matrices
    unity.fill(0), sigma_x.fill(0), sigma_y.fill(0), sigma_z.fill(0);
    mat dummy_unity={{1,0},{0,1}}, dummy_sigma_x={{0,1},{1,0}};
    mat dummy_sigma_y={{0,-1},{1,0}}, dummy_sigma_z={{1,0},{0,-1}};
    unity.set_real(dummy_unity);
    sigma_x.set_real(dummy_sigma_x);
    sigma_y.set_imag(dummy_sigma_y);
    sigma_z.set_real(dummy_sigma_z);

    // Rhombohedral has the shape of the unit cell, so don't need to carve out
    if (myParticle.shape == "Rhombohedral") {
        sizeCells = pSizeCells; 
    } 

    // Set up lattice and calculate number of cells within particle
    for (int i=-sizeCells; i<sizeCells; i++) {
        for (int j=-sizeCells;  j<sizeCells; j++) {
            for (int k=-sizeCells; k<sizeCells; k++) {
                dummyVec=i*latVec1+j*latVec2+k*latVec3-comShift;
                if (myParticle.WithinParticle(dummyVec)) {
                    pcomShift=pcomShift+2*dummyVec+sublatVec;
                    pCells++;
                }
            }
        }
    }   

    // Finding particle's COM
    pcomShift=0.5*pcomShift/pCells;

    // Set up data structures for the particle
    myParticle.SetDataStructures(pCells);

    // Carve out the unit sphere out of lattice
    for (int i=-sizeCells; i<sizeCells; i++) {
        for (int j=-sizeCells;  j<sizeCells; j++) {
            for (int k=-sizeCells; k<sizeCells; k++) {
                dummyVec=i*latVec1+j*latVec2+k*latVec3-comShift;
                if (myParticle.WithinParticle(dummyVec)) {
                    myParticle.sublatOne(index,0)=index;
                    myParticle.sublatOne(index,span(1,3))=dummyVec-pcomShift; 
                    myParticle.sublatTwo(index,0)=index;
                    myParticle.sublatTwo(index,span(1,3))=dummyVec+sublatVec-pcomShift;
                    index++;
                    if (myParticle.OutsideOrderCore(dummyVec)) {
                        dummy++;
                    }
                }
            }
        }
    }   

    dummy=0; index=0;
    // Set up the Hamiltonian
    for (int i=0; i<pCells; i++) {
        for (int j=0; j<pCells; j++) {

            /* dummyVec=sublatOne(i,span(1,3))-sublatTwo(j,span(1,3));*/
            dummyVec=myParticle.sublatOne(i,span(1,3))-myParticle.sublatTwo(j,span(1,3));

            // Set up nearest neighbour hoppings (between sublattices)
            if (dot(dummyVec,dummyVec) < sublatVecNorm*sublatVecNorm+delta) {
                if (i==j) {
                    // hoppings within a cell i
                    myParticle.Hamiltonian(span(4*i+0,4*i+1),span(4*i+2,4*i+3))+=(tHopping+deltatHopping)*unity;
                    myParticle.Hamiltonian(span(4*i+2,4*i+3),span(4*i+0,4*i+1))+=(tHopping+deltatHopping)*unity;
                } else {
                    // hoppings from sublatTwo of cell j to sublatOne of cell i
                    myParticle.Hamiltonian(span(4*i+0,4*i+1),span(4*j+2,4*j+3))+=(tHopping)*unity;

                    // hoppings from sublatOne of cell i to sublatTwo of cell j
                    myParticle.Hamiltonian(span(4*j+2,4*j+3),span(4*i+0,4*i+1))+=(tHopping)*unity;
                }
            }

            // Set up next nearest neighbour hoppings (within each sublattice)
            dummyVec=myParticle.sublatOne(i,span(1,3))-myParticle.sublatOne(j,span(1,3)); 
            if ( (dot(dummyVec,dummyVec) < latVecNorm*latVecNorm+delta)  && 
                    (dot(dummyVec,dummyVec) > sublatVecNorm*sublatVecNorm+delta) ) { 

                // Matrix elements are direction-dependent
                if ( (dummyVec(0) > 0 && dummyVec(1) > 0) || (dummyVec(0) > 0 &&
                            dummyVec(2) > 0) || (dummyVec(1) > 0 && dummyVec(2) > 0) ) {
                    dummyVec = cross(sublatVec, dummyVec-sublatVec);
                } else if ( (dummyVec(1) < 0 && dummyVec(2) < 0) || (dummyVec(0) > 0 &&
                            dummyVec(1) < 0) || (dummyVec(0) > 0 && dummyVec(2) < 0) ) {
                    dummyVec = cross(-latVec1+sublatVec, dummyVec+latVec1-sublatVec);
                } else if ( (dummyVec(0) < 0 && dummyVec(1) > 0) || (dummyVec(0) < 0 &&
                            dummyVec(2) < 0) || (dummyVec(1) > 0 && dummyVec(2) < 0) ) {
                    dummyVec = cross(-latVec2+sublatVec, dummyVec+latVec2-sublatVec);
                } else if ( (dummyVec(0) < 0 && dummyVec(2) > 0) || (dummyVec(1) < 0 &&
                            dummyVec(2) > 0) || (dummyVec(0) < 0 && dummyVec(1) < 0) ) {
                    dummyVec = cross(-latVec3+sublatVec, dummyVec+latVec3-sublatVec);
                } else {}
                dummyMat=dummyVec(0)*sigma_x+dummyVec(1)*sigma_y+dummyVec(2)*sigma_z;

                // See if both atoms are connected to the same linking atom
                // hoppings from cell j to cell i within sublattice One
                myParticle.Hamiltonian(span(4*i,4*i+1),span(4*j,4*j+1))=spinOrbitCoupling*dummyMat;

                // hoppings from cell i to cell j within sublattice One
                // Minus sign for hopping in opposite direction
                myParticle.Hamiltonian(span(4*j,4*j+1),span(4*i,4*i+1))=-spinOrbitCoupling*dummyMat;
                // hoppings from cell j to cell i within sublattice Two
                // Minus sign for sublattice  
                myParticle.Hamiltonian(span(4*i+2,4*i+3),span(4*j+2,4*j+3))=-spinOrbitCoupling*dummyMat;

                // hoppings from cell i to cell j within sublattice Two
                myParticle.Hamiltonian(span(4*j+2,4*j+3),span(4*i+2,4*i+3))=spinOrbitCoupling*dummyMat;

            }
        } 
    }

    //  myParticle.AddDisorder(disorderStrength, disorderCoverage, delta);

    // abs returns a matrix, max returns a vector, 2nd max - largest value
    if (max(max(abs(myParticle.Hamiltonian-myParticle.Hamiltonian.t()))) > delta ) {
        cout << "Error, Hamiltonian is not Hermitian by at least " << delta << endl;
    }

    // Get the eigenvectors and eigevalues
    eig_sym(myParticle.eigvals, myParticle.eigvecs, myParticle.Hamiltonian);

    myParticle.PrintInfo(disorderStrength);

    // Find the middle of the spectrum (Dirac point, E=0, for clean spectrum)
    for (int i=0; i<4*pCells-1;i++) {
        if (myParticle.eigvals(i) > 0 && myParticle.eigvals(i-1) < 0) {
            dummy/*flag*/ =i ;
        }
    }
    cout << "State above E #: " << dummy << " with E= " << myParticle.eigvals(dummy) << endl; 
    dummy=2*pCells;
    cout << "Picking state #: " << dummy << " with E= " << myParticle.eigvals(dummy) << endl; 

    mat probDensity(2*pCells,4+pCells); /* */
    probDensity(span(0,pCells-1),span(0,2))=myParticle.sublatOne.cols(1,3);
    probDensity(span(pCells,2*pCells-1),span(0,2))=myParticle.sublatTwo.cols(1,3);

    // Calculate probability density of probDensity labelled dummy and all probDensity above it at each lattice site
    for (int i=0; i<pCells; i++) {
        for (int j=0; j<1; j++) { /* */
            probDensity(i,3+j)=norm(myParticle.eigvecs(span(4*i,4*i+1),dummy+j))*norm(myParticle.eigvecs(span(4*i,4*i+1),dummy+j));
            probDensity(pCells+i,3+j)=norm(myParticle.eigvecs(span(4*i+2,4*i+3),dummy+j))*norm(myParticle.eigvecs(span(4*i+2,4*i+3),dummy+j));
        }    
    }   
    myParticle.eigvals.save("output_eigvals.txt",raw_ascii);
    probDensity.save("output_state.txt", raw_ascii);

    return 0;
}
