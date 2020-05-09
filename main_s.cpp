#include "declarations.h"

/*this function is simulating hydrodynamic interaction using the rotne
prager tensor and not the simple mobility function [X,Y]=full_hydro(x,
y, eta, a, dt, simtime, Kr, Kt, eps, rad, path) x,y are the starting
position nof the particles in the x, y axes.  eta is the viscosity, a
is the radius and dt is the time between collisions, simtime is the
number of steps for the simulation.  the output is x and y matrixs for
the positions */

// int main(int argc, char *argv[])
// {
//   int NParticles, MaxIter;
//   double Dist, aa, dt, K, Kt, Err, Eps, T;

//   NParticles = atoi(argv[1]); // Number of particles
//   Dist        = atof(argv[2]); // lattice parameter
//   Err        = atof(argv[3]); // dispersion in initail coords.
//   aa         = atof(argv[4]); // Particle diameter
//   dt         = atof(argv[5]); // Time step of simulation
//   MaxIter    = atoi(argv[6]); // Simulation length
//   K          = atof(argv[7]); // Radial confinment coeff.
//   Kt         = atof(argv[8]); // Tangential driving coeff.
//   Eps        = atof(argv[9]); // Wall coeff.
//   T          = atof(argv[10]);// Normalized temperature.
//   T=T*0.1;
//   // printf("T %f\n",T);
//   int ii, jj, NPrint;
//   double *D, *C, *R, *F, *P;
//   Point *PosC;
//   FILE* TrajFile, *CoordFile, *LogFile;
//   char filename[50],str[3];

//   Randomize();
//   long idum = -random();

//   //for(ii=0 ; ii < 100000 ; ii++) printf("%.16f\n",ran_nrc(&idum));
//   //for(ii=0 ; ii < 100000 ; ii++) printf("%.16f\n",Randn(idum));

//   /*#pragma omp parallel
//   {
// #pragma omp for
//     for(ii=0 ; ii < 10000 ; ii++) printf("%.16f\n",Randn(idum));
//     //for(ii=0 ; ii < 100000 ; ii++) printf("%d\n",rand());
//   }
//   exit(0);*/

//   InitSim(NParticles, Dist, Err, idum);

//   sprintf(str,"Traj%2d-%.1f-%.1f",NParticles,Kt,K);


//   LogFile = fopen("Log.txt","w");
//   TrajFile = fopen(str,"w");
// }
