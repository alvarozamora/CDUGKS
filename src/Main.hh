#ifndef MAIN_HH
#define MAIN_HH
#include <assert.h>
#include "Mesh.hh"

//Dimension and Resolution
int D = 3;
int N[3] = {128, 1, 1};       // For Lower D problem, set size to 1.
int NV[3] = {128, 1, 1};     //Both N and NV must be multiple of 4 (or 1 for lower D) -- (Newton-Cotes)
int Nc = N[0]*N[1]*N[2];    // Cells
int Nv = NV[0]*NV[1]*NV[2]; // Velocities
int effD = 1.0; //TODO: POTENTIAL BUG

//Physical Constants
double R = 0.5;         // Gas Constant
double K = 2.0;          //Internal DOF
double Cv = (3+K)*R/2;   //Specific Heat
double gma = (K+5)/(K+3); //gamma -- variable name taken
double w = 0.5;  //Viscosity exponent
double ur = 1.0; //Reference Visc
double Tr = 1.0; //Reference Temp

double Vmin[3] = {-10,0,0};
double Vmax[3] = {10,0,0};


void Cotes(double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ)
{
  int k,n;
  double dh, b, a;

  
  assert(NV[1]%4 == 0 || NV[1] == 1 && NV[2]%4 == 0 || NV[2] == 1); //Using np = 4 Newton Cotes


  //First Dimension
  a=Vmin[0]; b=Vmax[0];
  n=(NV[0])/4;
  dh=(b-a)/n; //old from dugks
  dh=(b-a)/(NV[0]-1)*4; //mine (Al)

  for(int kx = 0; kx < NV[0];kx++){
  	Co_X[kx]=Vmin[0]+kx*dh/4;
  } 
  for(int kx = 0;kx < n;kx++)
  {
    Co_WX[4*kx]=14.0;
    Co_WX[4*kx+1]=32.0;
    Co_WX[4*kx+2]=12.0;
    Co_WX[4*kx+3]=32.0;
  }
  Co_WX[0]=7.0;
  Co_WX[NV[0]-1]=7.0;

  //Check
  /*
  for(int kx = 0; kx < NV[0]; kx++){
    printf("Main.hh Co_X[%d] = %f\n", kx, Co_X[kx]);

  }
  for(int kx = 0; kx < NV[0]; kx ++){
    printf("Main.hh Co_WX[%d] = %f\n", kx, Co_WX[kx]);
  }
   for(int kx = 0; kx < NV[0]; kx++){Co_WX[kx]*=dh/90.;}
  */

  //Second Dimension
  a=Vmin[1]; b=Vmax[1];
  n=(NV[1]-1)/4;
  dh=(b-a)/n; //old from dugks
  dh=(b-a)/(NV[1]-1)*4; //mine (Al)

  for(int ky = 0; ky < NV[1]; ky++){
  	Co_Y[ky]=Vmin[1]+ky*dh/4;
  }

  if(NV[1] >= 4){
  for(int ky = 0; ky < n; ky++)
  {
    Co_WY[4*ky]=14.0;
    Co_WY[4*ky+1]=32.0;
    Co_WY[4*ky+2]=12.0;
    Co_WY[4*ky+3]=32.0;
  }
  Co_WY[0]=7.0;
  Co_WY[NV[1]-1]=7.0;

  for(int ky = 0; ky < NV[1]; ky++){Co_WY[ky]*=dh/90;}
  }else if(NV[1] == 1 && Vmin[1] == 0){Co_WY[0] = 1; Co_Y[0] = 0;}
  else{printf("Not Supported: NV[1]!= 1 && NV[1] < 4\n");}


  if(NV[2] >= 4){
  a=Vmin[2]; b=Vmax[2];
  n=NV[2]/4;
  dh=(b-a)/n; //old from dugks
  dh=(b-a)/(NV[2]-1)*4; //mine (Al)

  for(int kz = 0; kz < NV[2];kz++){
  	Co_Z[kz]=Vmin[2]+kz*dh/4;
  }
  for(int kz = 0; kz < n;kz++)
  {
    Co_WZ[4*kz]=14.0;
    Co_WZ[4*kz+1]=32.0;
    Co_WZ[4*kz+2]=12.0;
    Co_WZ[4*kz+3]=32.0;
  }
  Co_WZ[0]=7.0;
  Co_WZ[NV[2]-1]=7.0;

  for(int kz = 0; kz < N[2]; kz++){Co_WZ[kz]*=dh/90;}
  } else if(NV[2] == 1 && Vmin[2] == 0){Co_WZ[0] = 1; Co_Z[0] = 0;}
  else{printf("Not Supported: NV[2]!= 1 && NV[2] < 4\n");}
}


#endif