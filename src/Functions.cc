#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Functions.hh"

extern double* Co_X;
extern double* Co_WX;
extern double* Co_Y;
extern double* Co_WY;
extern double* Co_Z;
extern double* Co_WZ;


double Temperature(double E, double u){

	double T = (gma - 1)/R * (E - 0.5*u*u);

	return T;
}

double geq(int kx, int ky, int kz, double rho, double* U, double T, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ)
{
  double c2, x;
  double PI = 4.0*atan(1.0);
  c2 = sqrt( (Co_X[kx]-U[0])*(Co_X[kx]-U[0]) + (Co_Y[ky]-U[1])*(Co_Y[ky]-U[1]) + (Co_Z[kz]-U[2])*(Co_Z[kz]-U[2]));
  x  = Co_WX[kx]*Co_WY[ky]*Co_WZ[kx]*rho*exp(-c2/(2*R*T))/sqrt(2*PI*R*T)/sqrt(2*PI*R*T)/sqrt(2*PI*R*T);
  return x;
}