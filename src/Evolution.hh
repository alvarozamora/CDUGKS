#ifndef EVOLUTION_HH
#define EVOLUTION_HH
#include <assert.h>

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
extern double K;


double TimeStep();

void Evolve(double* g, double* b, double* gbar, double* bbar, double* gbarp, double* bbarp, double* Sg, double* Sb, double* rho, double* rhov, double* rhoE, int effD, double dt, double Tf, double Tsim, double dtdump, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double* gsigma, double* bsigma, double* gsigma2, double* bsigma2);

void Step1a(double* g, double* b, double* gbar, double* bbar, double* gbarp, double* bbarp, double* Sg, double* Sb, double* rho, double* rhov, double* rhoE, int effD, double dt, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ);
void Step1b(double* gbarp, double* bbarp, int effD, double* gsigma, double* bsigma, double* gsigma2, double* bsigma2);
void Step1c();

void Step2a();
void Step2b();
void Step2c();

void Step3();
void Step4();
void Step5();

#endif