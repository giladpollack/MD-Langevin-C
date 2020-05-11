#include <iostream>
#include <time.h>
#include <omp.h>
#include <assert.h>
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
  long InitVal = (long)-329475;
  InitVal -= time(0); // Randomizing the seed using the current time

  RandGen rng(&InitVal);
  rng.Randomize(); // Randomizing a bit more by calling the random function several times depending on time

  const int NumOfSims = 10;
  char Foldernames[NumOfSims][10];
  char BaseFoldername[] = "ParCheck";
  long Seeds[MAX_PARALLEL_SIMS] = {743023948,30974320,58624534,95832034,43709834,4309323,83470921,45894323,348749189,2389123};
  #pragma omp parallel for
  for (int i = 0; i < NumOfSims; i++)
  {
    sprintf(Foldernames[i],"%s %d",BaseFoldername, i + 1);
    Seeds[i] += time(0);
    RandGen currGen(&(Seeds[i]));
    InfoChamber(1e6, 1e-4, 20, 1e-6, 300, 0.001, 11e-6, 11e-6, 0, 7, Foldernames[i], false, currGen);
  }

  // Non-Parallel run
  // char SaveFoldername[20] = "Check";
  // InfoChamber(1e5, 1e-4, 20, 1e-6, 300, 0.001, 11e-6, 11e-6, 0, 7, SaveFoldername, false, rng);
}
