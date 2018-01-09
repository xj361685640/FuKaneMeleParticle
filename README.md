# Tight-binding model for topological nanoparticles

Microscopic tight-binding model applied to topological insulator nanoparticles as in:

G. Siroki, P. D. Haynes, D. K. K. Lee and V. Giannini Phys. Rev. Mater. 1, 024201 (2017)

The original tight-binding Hamiltonian is presented here:

L. Fu, C. L. Kane, and E. J. Mele, Phys. Rev. Lett. 98, 106803 (2007)

Spherical, Cubic, Rhombohedral and Spherical Shell nanoparticles are supported. These inherit common properties 
from the Particle class. An example input.txt file is provided. To compile run 'make'. 

## Getting Started
To compile type
```
make
```
To run execute
``` 
./main
```
which will read the file input.txt (example provided) for particle's shape and size and run the calculation.
The output files contain the eigenvalues and probability density of the associated wavefunctions.
### Dependencies

Armadillo, LAPACK, BLAS

## Authors

* **Gleb Siroki**

## License

This project is licensed under the MIT License - see the LICENSE.md file for details
