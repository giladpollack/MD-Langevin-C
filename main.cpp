#include <iostream>
#include <time.h>
#include <omp.h>
#include <assert.h>
#include <random>
#include <direct.h>
#include "declarations.h"
#include "Classes/RandGen.h"



const int MAX_PARALLEL_SIMS = 10;

/*this function is simulating hydrodynamic interaction using the rotne
prager tensor and not the simple mobility function [X,Y]=full_hydro(x,
y, eta, a, dt, simtime, Kr, Kt, eps, rad, path) x,y are the starting
position nof the particles in the x, y axes.  eta is the viscosity, a
is the radius and dt is the time between collisions, simtime is the
number of steps for the simulation.  the output is x and y matrixs for
the positions */
int main(int argc, char *argv[])
{  
  RandGen rng(time(NULL));

  const int NumOfSims = 1000;
  char Foldernames[NumOfSims][50];
  char BaseFoldername[] = "5 Particles With Interactions Shrunk 6 micron";
  int NumOfSimsRun = 0;
  
  // Creating the base folder to save the simulations in
  mkdir(BaseFoldername);

  // Running the parallel simulations
  #pragma omp parallel for
  for (int i = 0; i < NumOfSims; i++)
  {
    sprintf(Foldernames[i],"%s/%d",BaseFoldername, i + 1);
    RandGen currGen(rng.Randu(-1e7, 1e7));
    InfoChamber((long)1e7, 1e-4, 20, 1e-6, 300, 0.001, 11e-6, 11e-6, 6e-6, 5, Foldernames[i], false, currGen);
    NumOfSimsRun++;
    std::cout << "Run " << NumOfSimsRun << " out of " << NumOfSims << std::endl;
  }

  // Non-Parallel run
  // char SaveFoldername[20] = "Single system rand";
  // InfoChamber(1e8, 1e-5, 20, 1e-6, 300, 0.001, 11e-6, 11e-6, 0, 1, SaveFoldername, false, rng);
}