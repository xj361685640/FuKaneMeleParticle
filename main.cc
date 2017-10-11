#include <iostream>
#include <armadillo>                                                                                 
#include<complex>                                                                                    
#include<cstdlib>
#include <ctime>
using namespace std;                                                                                 
using namespace arma;                                                                                

class Particle {
    public:
        Particle() {}
        ~Particle() {}   

        int cells;
        int disorderedSites;
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

class Sphere: public Particle {
    private:
        int index;
        double dummy;
    public:
        double radius;
        double orderRadius;
        explicit Sphere(double pRadius, double orderpRadius) {
            radius = pRadius; // radius in distance units
            orderRadius = orderpRadius; // radius above which disorder starts
        }

        void AddDisorder(double disorderStrength, double disorderCoverage, double delta){ 
            rowvec dummyVec(3);
            disorderedSites=0;
            srand(time(0));//srand(3); 
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

};

int main(){
    // Input 
    int sizeCells;  //sizeCells of lattice in +/-ve x/y/z direction
    int pRadiusCells;  // pRadius of particle in unit cells - carved out of lattice
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
    while (inputfile >> sizeCells >> pRadiusCells >> disorderStrength >> disorderCoverage) 
        cout << "Input: "<< sizeCells << " " << pRadiusCells << " " <<  disorderStrength << " " << disorderCoverage << endl;

    Sphere mySphere(pRadiusCells*latVecNorm, pRadiusCells*latVecNorm - sublatVecNorm);

    int nCells = 8*sizeCells*sizeCells*sizeCells; // # unit cells 
    int pCells = 0; // # unit cells belonging to the particle

    int dummy = 0, index = 0; 
    rowvec dummyVec(3);
    cx_mat dummyMat(2,2);
    rowvec pcomShift(3);

    mat sublatOne(nCells, 4); // Sites at origin of a unit cell
    mat sublatTwo(nCells, 4); // Sites at sublatVec in a unit cell

    cx_mat unity(2,2), sigma_x(2,2), sigma_y(2,2), sigma_z(2,2); // Pauli matrices
    unity(0,0).real()=1; unity(1,1).real()=1;
    sigma_x(1,0).real()=1; sigma_x(0,1).real()=1;
    sigma_y(0,1).imag()=-1; sigma_y(1,0).imag()=1;
    sigma_z(0,0).real()=1; sigma_z(1,1).real()=-1;

    // Set up lattice and calculate number of cells within particle
    for (int i=-sizeCells; i<sizeCells; i++) {
        for (int j=-sizeCells;  j<sizeCells; j++) {
            for (int k=-sizeCells; k<sizeCells; k++) {
                dummyVec=i*latVec1+j*latVec2+k*latVec3-comShift;
                sublatOne(index,0)=index;
                sublatOne(index,span(1,3))=dummyVec;
                sublatTwo(index,0)=index;
                sublatTwo(index,span(1,3))=dummyVec+sublatVec;
                index++;
                if (dot(dummyVec,dummyVec) < mySphere.radius*mySphere.radius) {   
                    pcomShift=pcomShift+2*dummyVec+sublatVec;
                    pCells++;
                }
            }
        }
    }   

    // Finding particle's COM
    pcomShift=0.5*pcomShift/pCells;
    index=0;

    // Set up data structures for the particle
    mat copy_sublatOne(pCells, 4); // Used to visualise surface atoms - see the end
    mySphere.SetDataStructures(pCells);

    // Carve out the unit sphere out of lattice
    for (int i=-sizeCells; i<sizeCells; i++) {
        for (int j=-sizeCells;  j<sizeCells; j++) {
            for (int k=-sizeCells; k<sizeCells; k++) {
                dummyVec=i*latVec1+j*latVec2+k*latVec3-comShift;
                if (dot(dummyVec,dummyVec) < mySphere.radius*mySphere.radius) {  
                    copy_sublatOne(index,0)=index;
                    copy_sublatOne(index,span(1,3))=dummyVec;

                    mySphere.sublatOne(index,0)=index;
                    mySphere.sublatOne(index,span(1,3))=dummyVec-pcomShift; 
                    mySphere.sublatTwo(index,0)=index;
                    mySphere.sublatTwo(index,span(1,3))=dummyVec+sublatVec-pcomShift;
                    index++;
                    if (dot(dummyVec,dummyVec) > mySphere.orderRadius*mySphere.orderRadius) {
                        dummy++;
                    }
                }
            }
        }
    }   

    int pDisorderedCells=dummy;
    mat surf_state(2*pDisorderedCells,4);
    dummy=0; index=0;
    // Set up the Hamiltonian
    for (int i=0; i<pCells; i++) {
        for (int j=0; j<pCells; j++) {

            /* dummyVec=sublatOne(i,span(1,3))-sublatTwo(j,span(1,3));*/
            dummyVec=mySphere.sublatOne(i,span(1,3))-mySphere.sublatTwo(j,span(1,3));

            // Set up nearest neighbour hoppings (between sublattices)
            if (dot(dummyVec,dummyVec) < sublatVecNorm*sublatVecNorm+delta) {
                if (i==j) {
                    // hoppings within a cell i
                    mySphere.Hamiltonian(span(4*i+0,4*i+1),span(4*i+2,4*i+3))+=(tHopping+deltatHopping)*unity;
                    mySphere.Hamiltonian(span(4*i+2,4*i+3),span(4*i+0,4*i+1))+=(tHopping+deltatHopping)*unity;
                } else {
                    // hoppings from sublatTwo of cell j to sublatOne of cell i
                    mySphere.Hamiltonian(span(4*i+0,4*i+1),span(4*j+2,4*j+3))+=(tHopping)*unity;

                    // hoppings from sublatOne of cell i to sublatTwo of cell j
                    mySphere.Hamiltonian(span(4*j+2,4*j+3),span(4*i+0,4*i+1))+=(tHopping)*unity;
                }
            }

            // Set up next nearest neighbour hoppings (within each sublattice)
            dummyVec=mySphere.sublatOne(i,span(1,3))-mySphere.sublatOne(j,span(1,3)); 
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
                mySphere.Hamiltonian(span(4*i,4*i+1),span(4*j,4*j+1))=spinOrbitCoupling*dummyMat;

                // hoppings from cell i to cell j within sublattice One
                // Minus sign for hopping in opposite direction
                mySphere.Hamiltonian(span(4*j,4*j+1),span(4*i,4*i+1))=-spinOrbitCoupling*dummyMat;
                // hoppings from cell j to cell i within sublattice Two
                // Minus sign for sublattice  
                mySphere.Hamiltonian(span(4*i+2,4*i+3),span(4*j+2,4*j+3))=-spinOrbitCoupling*dummyMat;

                // hoppings from cell i to cell j within sublattice Two
                mySphere.Hamiltonian(span(4*j+2,4*j+3),span(4*i+2,4*i+3))=spinOrbitCoupling*dummyMat;

            }
        } 
    }

    mySphere.AddDisorder(disorderStrength, disorderCoverage, delta);

    // abs returns a matrix, max returns a vector, 2nd max - largest value
    if (max(max(abs(mySphere.Hamiltonian-mySphere.Hamiltonian.t()))) > delta ) {
        cout << "Error, Hamiltonian is not Hermitian by at least " << delta << endl;
    }

    cout << "Total # of sublattice sites: " << 2*mySphere.cells << " within radius " << mySphere.radius << endl;
    cout << "# of surface sublattice sites: " << 2*pDisorderedCells << " above radius " << mySphere.orderRadius << endl;
    cout << "# of disordered sublattice sites: " << mySphere.disorderedSites << " Strength: " << disorderStrength << " above radius " << mySphere.orderRadius << endl;
    
    // Get the eigenvectors and eigevalues
    eig_sym(mySphere.eigvals, mySphere.eigvecs, mySphere.Hamiltonian);

    // Find the middle of the spectrum (Dirac point, E=0, for clean spectrum)
    for (int i=0; i<4*pCells-1;i++) {
        if (mySphere.eigvals(i) > 0 && mySphere.eigvals(i-1) < 0) {
            dummy/*flag*/ =i ;
        }
    }
    cout << "State above E #: " << dummy << " with E= " << mySphere.eigvals(dummy) << endl; 
    dummy=2*pCells;
    cout << "Picking state #: " << dummy << " with E= " << mySphere.eigvals(dummy) << endl; 

    mat probDensity(2*pCells,4+pCells); /* */
    probDensity(span(0,pCells-1),span(0,2))=mySphere.sublatOne.cols(1,3);
    probDensity(span(pCells,2*pCells-1),span(0,2))=mySphere.sublatTwo.cols(1,3);

    // Calculate probability density of probDensity labelled dummy and all probDensity above it at each lattice site
    for (int i=0; i<pCells; i++) {
        for (int j=0; j<1; j++) { /* */
            probDensity(i,3+j)=norm(mySphere.eigvecs(span(4*i,4*i+1),dummy+j))*norm(mySphere.eigvecs(span(4*i,4*i+1),dummy+j));
            probDensity(pCells+i,3+j)=norm(mySphere.eigvecs(span(4*i+2,4*i+3),dummy+j))*norm(mySphere.eigvecs(span(4*i+2,4*i+3),dummy+j));
        }    
    }   
    mySphere.eigvals.save("output_eigvals.txt",raw_ascii);
    probDensity.save("output_state.txt", raw_ascii);

return 0;
}
