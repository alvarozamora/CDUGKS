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

extern double Cv;
extern double R;
extern double gma;

void TestProblem(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, int testProblem);

void SodShock(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, double rhoL = 1, double rhoR = 0.125, double PL= 1, double PR = 0.1); // rhoL, rhoR, PL, PR
void KHI();
void RTI();

#endif

//, double* rho, double* rhov, double* rhoE