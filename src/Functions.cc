#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Functions.hh"

extern double gma;
extern double K;
extern double R;
extern double ur;
extern double w;
extern double Tr;

double Temperature(double E, double u){

	double T = (gma - 1)/R * (E - 0.5*u*u);

	return T;
}

double geq(double vx, double vy, double vz, double rho, double* U, double T, double wx, double wy, double wz)
{
  double c2, x;
  double PI = 4.0*atan(1.0);
  c2 = sqrt( (vx-U[0])*(vx-U[0]) + (vy-U[1])*(vy-U[1]) + (vz-U[2])*(vz-U[2]));
  x  = wx*wy*wz*rho*exp(-c2/(2*R*T))/sqrt(2*PI*R*T)/sqrt(2*PI*R*T)/sqrt(2*PI*R*T); //TODO: POTENTIAL BUG (Dimensions)
  return x;
}

double visc(double T)
{
	return ur*pow(T/Tr,w);
}

double VanLeer(double L, double C, double R, double xL, double xC, double xR)
{
	double s1 = (C - L)/(xC - xL);
	double s2 = (R - C)/(xR - xC);

	return (sgn(s1) + sgn(s1))*(abs(s1) * abs(s2))/(abs(s1) + abs(s2));
}