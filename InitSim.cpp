#include "declarations.h"
#include "RandGen.h"

/* this function calculate the initial coordinates and parameters for a line array of particles.*/

void InitSim(int NParticles, double Dist, double Err, long &idum)
{
  int ii;
  Point *PosC;
  FILE* CoordFile;

  CoordFile = fopen("Coords.txt","w");
  //printf("init \n");

  PosC = (Point*) malloc (NParticles*sizeof(Point));
  if(PosC==NULL) exit(10);

  for (ii=0; ii < NParticles; ii++){
    PosC[ii].x   =  (1.0*ii*Dist) + RandGen::Randn()*Err;
    PosC[ii].y   =  RandGen::Randn()*Err;
  }

  // printf("initial positions \n");

  for(ii=0 ; ii < NParticles ; ii++) fprintf(CoordFile,"%.16f %.16f \n",PosC[ii].x,PosC[ii].y);

  fclose(CoordFile);
  free(PosC);
}
