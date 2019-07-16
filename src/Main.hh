#ifndef MAIN_HH
#define MAIN_HH
#include "Mesh.hh"

//Dimension and Resolution
int D = 3;
int N[3] = {128,1,1};       // For Lower D problem, set size to 1.
int NV[3] = {128, 1, 1};
int Nc = N[0]*N[1]*N[2];    // Cells
int Nv = NV[0]*NV[1]*NV[2]; // Velocities

//Physical Constants
double R = 0.5;
double K = 2.0;
double Cv = (3+K)*R/2;
double gma = (K+5)/(K+3); //gamma

#endif