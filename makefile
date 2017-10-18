CC=g++
CFLAGS=-Wall

main: main.cc Particle.cc Spherical.cc Cubic.cc Rhombohedral.cc
	$(CC) -std=c++11 $(CFLAGS) -o main main.cc Particle.cc Spherical.cc Shell.cc Cubic.cc Rhombohedral.cc -llapack -lblas -lgfortran -larmadillo


