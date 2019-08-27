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

double geq(double c2, double rho, double T)
{
  double x;
  double PI = 4.0*atan(1.0);
  //printf("vx = %f, vy = %f, vz %f, U0 = %f, U1 = %f, U2 = %f\n", vx, vy, vz, U[0], U[1], U[2]);

  x  = rho*exp(-c2/(2*R*T));
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
	if(xL > xC){ xL = xL - 1.0;}

	double s1 = (C - L)/(xC - xL);
	double s2 = (R - C)/(xR - xC);


	if(C - L == 0. && R - C == 0.){
		return 0;}

	else if(xL == xC || xR == xC){
		return 0;}
	else{
		//printf("L = %f, C = %f, R = %f, xL = %f, xL = %f, xR = %f, s1 = %f, s2 = %f, VL = %f\n", L, C, R, xL, xC, xR, s1, s2, (sgn(s1) + sgn(s2))*(abs(s1) * abs(s2))/(abs(s1) + abs(s2)));
		return (sgn(s1) + sgn(s2))*(abs(s1) * abs(s2))/(abs(s1) + abs(s2));}
}