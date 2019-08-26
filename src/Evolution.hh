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
extern double Pr;


void Evolve(double* g, double* b, double* gbar, double* bbar, double* gbarp, double* bbarp, double* Sg, double* Sb, double* rho, double* rhov, double* rhoE, double* dt, double Tf, double Tsim, double dtdump, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double* gsigma, double* bsigma, double* gsigma2, double* bsigma2, Cell* mesh, double* gbarpbound, double* bbarpbound, double* rhoh, double* rhovh, double* rhoEh);

void Step1a(double* g, double* b, double* gbar, double* bbar, double* gbarp, double* bbarp, double* Sg, double* Sb, double* rho, double* rhov, double* rhoE, double dt, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ);
void Step1b(double* gbarp, double* bbarp, double* gsigma, double* bsigma, double* gsigma2, double* bsigma2, Cell* mesh, double* gbarpbound, double* bbarpbound, double* Co_X, double* Co_Y, double* Co_Z);
void Step1c(double* gbar, double* bbar, double* gbarpbound, double*bbarpbound, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ, double* gsigma2, double* bsigma2, double dt);

void Step2a(double* gbar, double* bbar, double* Co_X, double* Co_Y, double* Co_Z, double* Co_WX, double* Co_WY, double* Co_WZ, double dt, double* rhoh, double* rhovh, double* rhoEh);
void Step2b(double* gbar, double* bbar, double dt, double* rhoh, double* rhovh, double* rhoEh, double* Co_X, double* Co_Y, double* Co_Z);
void Step2c(double* gbar, double* bbar, double* Co_X, double* Co_Y, double* Co_Z, Cell* mesh, double* Fg, double* Fb);

void Step3();
void Step4and5(double* rho, double* rhov, double* rhoE, double dt, Cell* mesh, double* Fg, double* Fb, double* Co_X, double* Co_Y, double* Co_Z, double* Co_WX, double* Co_WY, double* Co_WZ, double* g, double* b);

double TimeStep(double dt, double dtdump, double tend);

#endif