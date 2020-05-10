#include "declarations.h"
#include "RandGen.h"
#include <iostream>


/*this function is simulating hydrodynamic interaction using the rotne
prager tensor and not the simple mobility function [X,Y]=full_hydro(x,
y, eta, a, dt, simtime, Kr, Kt, eps, rad, path) x,y are the starting
position nof the particles in the x, y axes.  eta is the viscosity, a
is the radius and dt is the time between collisions, simtime is the
number of steps for the simulation.  the output is x and y matrixs for
the positions */
int main(int argc, char *argv[])
{
  // TODO: Fix to make the random better (use current time as seed)!
  long initVal = (long)-329475;
  RandGen::Init(&initVal);
  
  int NParticles, MaxIter;
  double Dist, aa, dt, K, Kt, Err, Eps, T;
  int numArgs = argc;

  // Basic definitions
  char SaveFoldername[20] = "Check";
  InfoChamber(1e6, 1e-4, 20, 1e-6, 300, 0.001, 11e-6, 11e-6, 0, 7, SaveFoldername, false);



  NParticles = atoi(argv[1]); // Number of particles
  Dist       = atof(argv[2]); // Vortex radius
  Err        = atof(argv[3]); // dispersion in initail coords.
  aa         = atof(argv[4]); // Particle diameter
  dt         = atof(argv[5]); // Time step of simulation
  MaxIter    = atoi(argv[6]); // Simulation length
  K          = atof(argv[7]); // Radial confinment coeff.
  Kt         = atof(argv[8]); // Tangential driving coeff.
  Eps        = atof(argv[9]); // Wall coeff.
  T          = atof(argv[10]);// Normalized temperature.
  T=T*0.1;
  // printf("T %f\n",T);
  int ii, jj, NPrint;
  double *D, *C, *R, *F, *P;
  Point *PosC;
  Point *PosP;
  FILE* TrajFile, *CoordFile, *LogFile;
  char filename[50],str[3];



  //for(ii=0 ; ii < 100000 ; ii++) printf("%.16f\n",ran_nrc(&idum));
  //for(ii=0 ; ii < 100000 ; ii++) printf("%.16f\n",Randn(idum));

  /*#pragma omp parallel
  {
#pragma omp for
    for(ii=0 ; ii < 10000 ; ii++) printf("%.16f\n",Randn(idum));
    //for(ii=0 ; ii < 100000 ; ii++) printf("%d\n",rand());
  }
  exit(0);*/
  // InitSim(NParticles, Dist, Err, idum);
  sprintf(filename,"Traj%2d-%.1f-%.1f",NParticles,Kt,K);
  LogFile = fopen("Log.txt","w");
  TrajFile = fopen(filename,"w");

  NPrint = 1000;

  /*****************************************/
  /* Read Coordinates and input parameters */
  /*****************************************/

  CoordFile = fopen("Coords.txt","r");

  PosC = (Point*) malloc (NParticles*sizeof(Point));
  if(PosC==NULL) exit(0);

  for(ii=0 ; ii < NParticles ; ii++) fscanf(CoordFile,"%lf %lf \n",&PosC[ii].x,&PosC[ii].y);

  fclose(CoordFile);

  /**********************/
  /* Memory allocations */
  /**********************/

  D = (double*) malloc (4*NParticles*NParticles*sizeof(double));
  if(D==NULL) exit(0);
  C = (double*) malloc (4*NParticles*NParticles*sizeof(double));
  if(C==NULL) exit(1);
  R = (double*) malloc (2*NParticles*sizeof(double));
  if(R==NULL) exit(2);
  F = (double*) malloc (2*NParticles*sizeof(double));
  if(F==NULL) exit(3);
  P = (double*) malloc (2*NParticles*sizeof(double));
  if(P==NULL) exit(4);
  PosP = (Point*) malloc (NParticles*sizeof(Point));
  if(PosP==NULL) exit(5);

  fprintf(LogFile,"Started Simulation\n");
  fflush(LogFile);

  for(ii=2; ii < MaxIter ; ii++){

    // HydroRot(PosC,PosP,D,C,R,F,aa,dt,K,Kt,Eps,P,NParticles,T,idum);

    if(!(ii%NPrint)){

      for(jj=0 ; jj < NParticles ; jj++) fprintf(TrajFile,"%.16f %.16f %d\n",PosC[jj].x,PosC[jj].y,ii/NPrint);

      fflush(TrajFile);
    }
  }
  
  fprintf(LogFile,"Finished Simulation\n");
  fflush(LogFile);

  fclose(LogFile);
  fclose(TrajFile);
  //printf("closed file\n");

  free(PosC);
  free(D);
  free(C);
  free(R);
  free(F);
  free(P);
  free(PosP);
}
