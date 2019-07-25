#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH


double Temperature(double E, double u);

//double geq(int kx, int ky, int kz, double rho, double* U, double T, double* Co_X, double* Co_WX, double* Co_Y, double* Co_WY, double* Co_Z, double* Co_WZ);

double geq(double vx, double vy, double vz, double rho, double* U, double T, double wx, double wy, double wz);

#endif