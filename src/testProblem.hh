#ifndef TESTPROBLEM_HH
#define TESTPROBLEM_HH

#include "Mesh.hh"
#include "Functions.hh"

extern int testProblem;
extern int D;
extern int Nc;
extern int Nv;
extern int N[3];
extern int NV[3];
extern int effD;

extern double Cv;
extern double R;
extern double gma;
extern double K;


void TestProblem(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, int testProblem, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ);

void SodShock(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double rhoL = 1, double rhoR = 0.125, double PL= 1, double PR = 0.1); // rhoL, rhoR, PL, PR
void KHI();
void RTI();

#endif