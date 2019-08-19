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
extern int effD;

double Temperature(double E, double u){

	double T = (gma - 1)/R * (E - 0.5*u*u);

	return T;
}

double geq(double vx, double vy, double vz, double rho, double* U, double T, double wx, double wy, double wz)
{
  double c2, x;
  double PI = 4.0*atan(1.0);
  //c2 = sqrt( (vx-U[0])*(vx-U[0]) + (vy-U[1])*(vy-U[1]) + (vz-U[2])*(vz-U[2])); why was this square rooted?
  c2 = (vx-U[0])*(vx-U[0]) + (vy-U[1])*(vy-U[1]) + (vz-U[2])*(vz-U[2]);
  //printf("vx = %f, vy = %f, vz %f, U0 = %f, U1 = %f, U2 = %f\n", vx, vy, vz, U[0], U[1], U[2]);
  x  = wx*wy*wz*rho*exp(-c2/(2*R*T)); //TODO: POTENTIAL BUG (Dimensions)
  for(int dim = 0; dim < effD; dim++){x = x/sqrt(2*PI*R*T);}
  return x;
}

double visc(double T)
{
	return ur*pow(T/Tr,w);
}

double VanLeer(double L, double C, double R, double xL, double xC, double xR)
{
	if(xR < xC){ xR = xR + 1.0;}
	if(xL > xR){ xL = xL - 1.0;}

	double s1 = (C - L)/(xC - xL);
	double s2 = (R - C)/(xR - xC);


	if(C - L == 0. && R - C == 0.)
		{return 0;}
	else{
		//printf("L = %f, C = %f, R = %f, xL = %f, xL = %f, xR = %f, s1 = %f, s2 = %f, VL = %f\n", L, C, R, xL, xC, xR, s1, s2, (sgn(s1) + sgn(s2))*(abs(s1) * abs(s2))/(abs(s1) + abs(s2)));
		return (sgn(s1) + sgn(s2))*(abs(s1) * abs(s2))/(abs(s1) + abs(s2));}
}