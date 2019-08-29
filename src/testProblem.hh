#ifndef TESTPROBLEM_HH
#define TESTPROBLEM_HH

#include "Mesh.hh"
#include "Functions.hh"

/*
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
*/

void TestProblem(int* N, int* NV, int* Nc, int* Nv, int* BCs, double* Vmin, double* Vmax, int testProblem, double* R, double* K, double* Cv, double* gma, double* w , double* ur, double* Tr, double* Pr, int* effD);
void InitializeTestProblem(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, int testProblem, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double R, double K, double Cv, double gma, double w, double ur, double Tr, double Pr, int* N, int* NV, int effD);

//void SodShock(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double rhoL = 1.0, double rhoR = 0.125, double PL= 1, double PR = 0.1); // rhoL, rhoR, PL, PR
void SodShock(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double R, double K, double Cv, double gma, double w, double ur, double Tr, double Pr, int* N, int* NV, int effD, double rhoL = 1.0, double rhoR = 0.125, double PL = 1, double PR = 0.1); // rhoL, rhoR, PL, PR
void KHI(Cell* mesh, double* g, double* b, double* rho, double* rhov, double* rhoE, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double R, double K, double Cv, double gma, double w, double ur, double Tr, double Pr, int* N, int* NV, int effD, double rhoT = 2.0, double rhoB = 1.0, double PT = 1.0, double PB = 1.0, double vrel = 2.0, double amp = 0.05);
void RTI();

#endif